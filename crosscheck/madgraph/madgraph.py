import subprocess
from particle import Particle
import warnings
import glob
import re
from smpl_io import io as sio
from smpl_io import sed

def pdgid_to_name(id):
    id = int(id)
    if id == 0:
        return 'g'
    if id == 22:
        return 'A'
    return Particle.from_pdgid(id).programmatic_name

def compute(in_df):
    global_diags = []
    procs = []
    # separate same flavour diagrams 
    # get every column that is prefixed with 'in_'
    in_cols = sorted([col for col in in_df.columns if str(col).startswith('in_')] )
    out_cols = sorted([col for col in in_df.columns if str(col).startswith('out_')])
    if len(in_cols) + len(out_cols) > 9: # TODO improve this
        warnings.warn('More than 9 particles, this is not supported')
    for i,row in in_df.iterrows():
        # get the ids of the particles
        in_ids = [row[col] for col in in_cols]
        out_ids = [row[col] for col in out_cols]
        # do not append to global_diags if it is already there
        if not (*in_ids, *out_ids) in global_diags:
            global_diags.append((*in_ids, *out_ids))
            cid = len(global_diags) - 1
            genstr = "import model sm\n"
            genstr += "generate" + " " + " ".join([pdgid_to_name(p) for p in in_ids]) + " > " + " ".join([pdgid_to_name(p) for p in out_ids]) + " [virt=QCD]\n"
            genstr += f"output /tmp/output_{cid}\n"
            genstr += "launch\n"
            print("Defining process: ", genstr, " with id: ", cid, " in madgraph.py")
            with open(f"/tmp/mg_{cid}", "w") as f:
                f.write(genstr)
            procs.append(subprocess.Popen(["/opt/original_MG/MG5_aMC_v3_5_0/bin/mg5_aMC", f"/tmp/mg_{cid}"]))
    # wait for all processes to finish
    for proc in procs:
        proc.wait()
    for cid in range(len(global_diags)):
        # replace READPS = .FALSE. with READPS = .TRUE.
        fn = glob.glob(f"/tmp/output_{cid}/SubProcesses/P*/check_sa.f")[0]
        sio.write(fn,sed.sed("READPS = .FALSE.", "READPS = .TRUE.",fn).read())
        fm = glob.glob(f"/tmp/output_{cid}/SubProcesses/P*/")[0]
        procs.append(subprocess.Popen(["make", "-C", fm, "check"]))
    for proc in procs:
        proc.wait()


    in_df['Born'] = 0.0
    in_df['Finite'] = 0.0
    in_df['Single Pole'] = 0.0
    in_df['Double Pole'] = 0.0

    for i,row in in_df.iterrows():
        in_ids = row[in_cols]
        out_ids = row[out_cols]
        cid = global_diags.index((*in_ids, *out_ids))
        pn_cols = sorted([col for col in in_df.columns if str(col).startswith('p_')] )
        psp = row[pn_cols].to_numpy().astype(float).reshape((int(len(pn_cols)/4),4)).tolist()
        with open(glob.glob(f"/tmp/output_{cid}/SubProcesses/P*/")[0] + "PS.input",'w') as f:
            for ps in psp:
                print(*ps, file=f)
        output = subprocess.check_output([glob.glob(f"/tmp/output_{cid}/SubProcesses/P*/check")[0]], cwd=glob.glob(f"/tmp/output_{cid}/SubProcesses/P*/")[0]).decode('utf-8')
        print(output)
        born = re.search(r"born\s*=\s*(.*)GeV", output).group(1)
        finite = re.search(r"finite\s*=\s*(.*)GeV", output).group(1)
        single = re.search(r"1eps\s*=\s*(.*)GeV", output).group(1)
        double = re.search(r"2eps\s*=\s*(.*)GeV", output).group(1)
        alphas = re.search(r"aS\s*=\s*(.*)", output).group(1)
        alpha1 = re.search(r"aEWM1\s*=\s*(.*)", output).group(1)
        in_df.loc[i,('Born',)] = float(born)
        in_df.loc[i,('Finite',)] = float(finite)
        in_df.loc[i,('Single Pole',)] = float(single)
        in_df.loc[i,('Double Pole',)] = float(double)
        in_df.loc[i,('alphas',)] = float(alphas)
        in_df.loc[i,('alpha',)] = 1/float(alpha1)
    return in_df
