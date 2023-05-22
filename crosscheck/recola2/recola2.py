from pyrecola import define_process_rcl,generate_processes_rcl,compute_process_rcl
from particle import Particle
import warnings

def pdgid_to_name(id):
    id = int(id)
    if id == 0:
        return 'g'
    if id == 22:
        return 'A'
    return Particle.from_pdgid(id).programmatic_name

def compute(in_df,order='NLO'):
    global_diags = []
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

            proc = (' '.join([pdgid_to_name(p) for p in in_ids]) 
                               + ' -> ' 
                               + ' '.join([pdgid_to_name(p) for p in out_ids]))
            print("Defining process: ", proc, " with id: ", cid, " and order: ", order, " in recola2.py")
            define_process_rcl(cid, proc , order)
    print("Generating processes in recola2.py")
    # generate it
    generate_processes_rcl()

    in_df['A0'] = 0.0
    in_df['A1'] = 0.0
    in_df['A2'] = 0.0
    for i,row in in_df.iterrows():
        in_ids = row[in_cols]
        out_ids = row[out_cols]
        cid = global_diags.index((*in_ids, *out_ids))
        pn_cols = sorted([col for col in in_df.columns if str(col).startswith('p_')] )
        # get the momenta of the particles
        psp = row[pn_cols].to_numpy().astype(float).reshape((int(len(pn_cols)/4),4)).tolist()
        #print(psp.to_numpy().astype(float))
        print("Computing process: ", cid, " with momenta: ", psp, " and order: ", order, " in recola2.py")
        A0, A1, A2 = compute_process_rcl(cid, psp, order)
        in_df.loc[i,('A0',)] = A0
        in_df.loc[i,('A1',)] = A1
        in_df.loc[i,('A2',)] = A2
    return in_df


