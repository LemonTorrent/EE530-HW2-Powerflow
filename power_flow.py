"""
Python version of the Chapman "power_flow" Matlab script.
Modified to support PQ, Constant Impedance (Z/LZ), and Constant Current (I/LI) load models.

Enter into Bash Command Line:
    python power_flow.py <input_filename>  (Use the given file that comes with the zip file package)

Input file format (same as Chapman):
  SYSTEM name baseMVA tol
  BUS name type volts Pgen Qgen Pload Qload Qcap
  LINE from to Rse Xse Gsh Bsh Rating

Supported bus types:
  SL = slack bus
  PV = generator bus with fixed voltage magnitude
  PQ / LP = constant-power load bus
  Z / LZ = constant-impedance load bus (P, Q vary with |V|^2)
  I / LI = constant-current load bus (P, Q vary with |V|)

Notes:
 - Per-unit arithmetic is used the same as the original Matlab script.
 - The solver uses the same update rules:
    * PV buses: Q updated each iteration and voltage magnitude kept fixed
    * SL (slack) bus: voltage angle fixed at 0 and magnitude fixed
 - Output mirrors the Matlab report: bus info, line flows, losses, alerts.

"""
import sys
import math
import cmath

def normalize_bus_type(btype):
    """Normalize supported bus-type aliases used in the input files."""
    btype = btype.upper()
    if btype in ('PQ', 'LP'):
        return 'PQ'
    if btype in ('Z', 'LZ'):
        return 'LZ'
    if btype in ('I', 'LI'):
        return 'LI'
    return btype

def read_system(OOEP_Initial):
    system = {}
    bus = []
    line = []
    with open(OOEP_Initial, 'r') as f:
        for raw in f:
            data = raw.strip()
            if data == '' or data.startswith('%'):
                continue
            tokens = data.split()
            key = tokens[0].upper()
            if key == 'SYSTEM':
                system['name'] = tokens[1]
                system['baseMVA'] = float(tokens[2])
                system['v_tol'] = float(tokens[3])
            elif key == 'BUS':
                # BUS name type volts Pgen Qgen Pload Qload Qcap
                name = tokens[1]
                btype = normalize_bus_type(tokens[2])
                # If Vbus given as two values (mag,angle) handle; otherwise single magnitude (angle=0)
                if len(tokens) >= 4:
                    Vbus_val = float(tokens[3])
                else:
                    Vbus_val = 1.0
                # Remaining numeric fields (some may be omitted; default to 0)
                nums = [float(x) for x in tokens[4:]]
                while len(nums) < 5:
                    nums.append(0.0)
                Pgen = nums[0] / system['baseMVA']
                Qgen = nums[1] / system['baseMVA']
                Pload = nums[2] / system['baseMVA']
                Qload = nums[3] / system['baseMVA']
                Qcap  = nums[4] / system['baseMVA']
                bus.append({
                    'name': name,
                    'type': btype,
                    'Vbus': Vbus_val,   # initial magnitude guess; solver will create complex values
                    'PG': Pgen,
                    'QG': Qgen,
                    'PL_nominal': Pload,
                    'QL_nominal': Qload,
                    'PL': Pload,
                    'QL': Qload,
                    'QC': Qcap,
                    'n_lines': 0
                })
            elif key == 'LINE':
                # LINE from to Rse Xse Gsh Bsh Rating
                from_name = tokens[1]
                to_name = tokens[2]
                nums = [float(x) for x in tokens[3:]]
                while len(nums) < 5:
                    nums.append(0.0)
                Rse, Xse, Gsh, Bsh, rating = nums[0], nums[1], nums[2], nums[3], nums[4]
                line.append({
                    'from_name': from_name,
                    'to_name': to_name,
                    'Rse': Rse,
                    'Xse': Xse,
                    'Gsh': Gsh,
                    'Bsh': Bsh,
                    'rating': rating,
                    'from': 0,
                    'to': 0
                })
            else:
                # ignore unknown lines but warn
                print(f"WARNING: Unknown card: {data}")

    # Link lines to bus indices and count bus connections
    n_bus = len(bus)
    for li in line:
        found_from = False
        found_to = False
        for idx, b in enumerate(bus):
            if li['from_name'] == b['name']:
                li['from'] = idx
                b['n_lines'] += 1
                found_from = True
            if li['to_name'] == b['name']:
                li['to'] = idx
                b['n_lines'] += 1
                found_to = True
        if not found_from:
            raise ValueError(f"Invalid from bus name: {li['from_name']}")
        if not found_to:
            raise ValueError(f"Invalid to bus name: {li['to_name']}")

    # Check for slack bus count
    sl_count = sum(1 for b in bus if b['type'] == 'SL')
    if sl_count == 0:
        raise ValueError("No slack bus specified!")
    if sl_count > 1:
        raise ValueError("More than one slack bus specified!")

    # Check generator buses types
    for idx, b in enumerate(bus):
        btype = normalize_bus_type(b['type'])
        if (b['PG'] != 0 or b['QG'] != 0) and (btype not in ('SL','PV')):
            raise ValueError(f"Generator bus {idx+1} specified as type {b['type']}")

    # i_line count approximate
    print("Input summary statistics:")
    print(f"{len(line) + len(bus) + 1:4d} lines in system file (approx)")
    print(f"{1:4d} SYSTEM lines")
    print(f"{len(bus):4d} BUS lines")
    print(f"{len(line):4d} LINE lines\n")

    return bus, line, system

def build_ybus(bus, line):
    n_bus = len(bus)
    ybus = [[0+0j for _ in range(n_bus)] for __ in range(n_bus)]
    for li in line:
        fr = li['from']
        to = li['to']
        Z = complex(li['Rse'], li['Xse'])
        if Z == 0:
            Yse = 0+0j
        else:
            Yse = 1.0 / Z
        Ysh = complex(li['Gsh'], li['Bsh'])
        ybus[fr][fr] += Yse + Ysh
        ybus[to][to] += Yse + Ysh
        ybus[fr][to] += -Yse
        ybus[to][fr] += -Yse
    # convert to python lists of complex (already)
    return ybus

def solve_system(bus, ybus, system):
    n_bus = len(bus)
    base = system['baseMVA']

    # Initialize complex voltage vector Vbus
    Vbus = [0+0j]*n_bus
    for i, b in enumerate(bus):
        # initial guess: real value = magnitude (angle 0)
        mag = b['Vbus']
        Vbus[i] = complex(mag, 0.0)

    acc_fac = 1.0
    eps = 1e-4
    n_iter = 0

    # Convert ybus to 2D list for indexing
    while True:
        n_iter += 1
        V_old = Vbus.copy()

        for i in range(n_bus):
            btype = normalize_bus_type(bus[i]['type'])
            if btype == 'SL':
                continue  # skip slack

            # Dynamic load modeling based on bus voltage magnitude
            v_mag = abs(Vbus[i])
            if btype == 'LZ':
                # Constant Impedance: P and Q proportional to |V|^2
                bus[i]['PL'] = bus[i]['PL_nominal'] * (v_mag ** 2)
                bus[i]['QL'] = bus[i]['QL_nominal'] * (v_mag ** 2)
            elif btype == 'LI':
                # Constant Current: P and Q proportional to |V|
                bus[i]['PL'] = bus[i]['PL_nominal'] * v_mag
                bus[i]['QL'] = bus[i]['QL_nominal'] * v_mag
            else:
                # Constant Power (PQ / LP or PV): P and Q constant
                bus[i]['PL'] = bus[i]['PL_nominal']
                bus[i]['QL'] = bus[i]['QL_nominal']

            # Update net P and Q injections
            bus[i]['P'] = bus[i]['PG'] - bus[i]['PL']
            if btype != 'PV':
                bus[i]['Q'] = bus[i]['QG'] - bus[i]['QL'] + bus[i]['QC']

            # If PV bus, update Q estimate from current voltages
            if btype == 'PV':
                temp = 0+0j
                for j in range(n_bus):
                    temp += ybus[i][j] * Vbus[j]
                temp = Vbus[i].conjugate() * temp
                bus[i]['Q'] = -temp.imag

            # Sum contributions from other buses
            temp = 0+0j
            for j in range(n_bus):
                if j != i:
                    temp -= ybus[i][j] * Vbus[j]

            # Add current injected term (P - jQ)/conj(V)
            # Note bus[i]['P'], ['Q'] are per-unit
            denom = Vbus[i].conjugate()
            if abs(denom) == 0:
                # small protection
                denom = 1e-9
            inj = (bus[i]['P'] - 1j*bus[i]['Q']) / denom
            temp = temp + inj

            Vnew = temp / ybus[i][i]
            # relaxation
            Vbus[i] = V_old[i] + acc_fac * (Vnew - V_old[i])

            # For PV bus enforce magnitude constant
            if btype == 'PV':
                mag_fixed = abs(V_old[i])
                if abs(Vbus[i]) != 0:
                    Vbus[i] = Vbus[i] * (mag_fixed / abs(Vbus[i]))

        # check convergence (both real and imag parts)
        diffs = [Vbus[k] - V_old[k] for k in range(n_bus)]
        maxdiff = 0.0
        for d in diffs:
            maxdiff = max(maxdiff, abs(d.real), abs(d.imag))
        if maxdiff < eps:
            break
        if n_iter > 5000:
            raise RuntimeError("Power flow did not converge after 5000 iterations")

    # Save back voltages
    for i in range(n_bus):
        bus[i]['Vbus'] = Vbus[i]

    return bus, n_iter

def r2p(c):
    """ rectangular complex to (mag, deg) """
    mag = abs(c)
    ang = math.degrees(cmath.phase(c))
    return mag, ang

def report_system(bus, line, system, ybus, n_iter):
    n_bus = len(bus)
    base = system['baseMVA']
    print("\nResults for Case", system.get('name','(unknown)'))
    print("|====================Bus Information====================|=======Line Information================|")
    print("Bus  Name       Type   Volts(pu) / angle(deg)   Gen(MW) Gen(MVAR)   Load(MW) Load(MVAR) Cap(MVAR)   ToBus    P(MW)   Q(MVAR)   S(MVA)")
    print("-"*120)
    PGtot = QGtot = PLtot = QLtot = QCtot = 0.0

    for i in range(n_bus):
        b = bus[i]
        mag, phase = r2p(b['Vbus'])
        if b['type'].upper() == 'SL':
            temp = 0+0j
            for j in range(n_bus):
                temp += ybus[i][j] * bus[j]['Vbus']
            temp = b['Vbus'].conjugate() * temp
            b['P'] = temp.real
            b['Q'] = -temp.imag

        Pgen = (b['P'] + b['PL']) * base
        Qgen = (b['Q'] + b['QL'] - b['QC']) * base
        PL = b['PL'] * base
        QL = b['QL'] * base
        QC = b['QC'] * base
        PGtot += Pgen
        QGtot += Qgen
        PLtot += PL
        QLtot += QL
        QCtot += QC

        print(f"{i+1:3d}  {b['name']:8s} {b['type']:3s}   {mag:6.3f}/ {phase:7.2f}    {Pgen:8.2f} {Qgen:9.2f}   {PL:8.2f} {QL:8.2f} {QC:8.2f}", end='')

        # list lines connected and flows
        count = 0
        for jj, li in enumerate(line):
            if li['from'] == i or li['to'] == i:
                # determine other end
                if li['from'] == i:
                    other = li['to']
                    other_name = li['to_name']
                else:
                    other = li['from']
                    other_name = li['from_name']

                # current on branch as in Matlab: il = (V_i - V_other) * ybus[i][other
                il = (b['Vbus'] - bus[other]['Vbus']) * ybus[i][other]
                pl = b['Vbus'] * il.conjugate()
                P_flow = -pl.real * base
                Q_flow = -pl.imag * base
                S_flow = math.sqrt(P_flow**2 + Q_flow**2)
                print(f"   {other_name:8s}   {P_flow:8.2f}  {Q_flow:8.2f}  {S_flow:8.2f}")
                count += 1
        if count == 0:
            print()

    print("-"*120)
    print(f"Totals: PG={PGtot:.2f} MW, QG={QGtot:.2f} MVAR, PL={PLtot:.2f} MW, QL={QLtot:.2f} MVAR, QC={QCtot:.2f} MVAR")
    print("-"*120)

    # Line losses
    print("\n                      Line Losses")
    print("| Line  From      To        Ploss(MW)  Qloss(MVAR) |")
    print("-"*55)
    Pltot = 0.0
    Qltot = 0.0
    for idx, li in enumerate(line):
        jj = li['from']
        kk = li['to']
        il = (bus[jj]['Vbus'] - bus[kk]['Vbus']) * ybus[jj][kk]
        Pl = (abs(il)**2) * li['Rse'] * base
        Ql = (abs(il)**2) * li['Xse'] * base
        Pltot += Pl
        Qltot += Ql
        print(f"{idx+1:5d}  {li['from_name']:8s}  {li['to_name']:8s}  {Pl:10.2f}  {Ql:10.2f}")
    print("-"*55)
    print(f"Totals: Ploss={Pltot:.2f} MW, Qloss={Qltot:.2f} MVAR")
    print("-"*55)

    # Alerts: voltages out of tolerance and line ratings
    print("\n                        Alerts")
    print("-"*55)
    n_alert = 0
    for i, b in enumerate(bus):
        if (abs(b['Vbus']) < (1.0 - system['v_tol'])) or (abs(b['Vbus']) > (1.0 + system['v_tol'])):
            n_alert += 1
            print(f"ALERT: Voltage on bus {b['name']} out of tolerance. |V|={abs(b['Vbus']):.4f}")

    for idx, li in enumerate(line):
        jj = li['from']
        kk = li['to']
        il = (bus[jj]['Vbus'] - bus[kk]['Vbus']) * ybus[jj][kk]
        Sl = abs(bus[jj]['Vbus'] * il.conjugate()) * base
        if Sl > li['rating']:
            n_alert += 1
            print(f"ALERT: Rating on line {idx+1} exceeded: {Sl:.2f} MVA > {li['rating']:.2f} MVA")

    if n_alert == 0:
        print("NONE")
    print(f"\nDone in {n_iter} iterations.\n")

def main(argv):
    if len(argv) < 2:
        print("Usage: python power_flow.py OOEP_Initial")
        return
    OOEP_Initial = argv[1]
    bus, line, system = read_system(OOEP_Initial)
    ybus = build_ybus(bus, line)
    bus, n_iter = solve_system(bus, ybus, system)
    report_system(bus, line, system, ybus, n_iter)

if __name__ == '__main__':
    main(sys.argv)
