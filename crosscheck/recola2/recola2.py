from pyrecola import *
from particle import Particle
from math import pi
import tqdm
import warnings
set_delta_ir_rcl(0,pi**2/6)
set_dynamic_settings_rcl(0)

# Let's print the squared amplitude
set_print_level_squared_amplitude_rcl(1)
set_dynamic_settings_rcl(1)


def pdgid_to_name(id):
    id = int(id)
    if id == 0:
        return 'g'
    if id == 22:
        return 'A'
    return Particle.from_pdgid(id).programmatic_name.replace("_bar","~")

def compute(in_df,order='NLO',debug=False):
    global_diags = []
    if debug:
        # The standard output is selected
        set_output_file_rcl('*')
    # separate same flavour diagrams 
    # get every column that is prefixed with 'in_'
    in_cols = sorted([col for col in in_df.columns if str(col).startswith('in_')] )
    out_cols = sorted([col for col in in_df.columns if str(col).startswith('out_')])
    if len(in_cols) + len(out_cols) > 9: # TODO improve this
        warnings.warn('More than 9 particles, this is not supported')
    for i,row in tqdm.tqdm(in_df.iterrows(),desc="Defining processes in recola2.py",total=in_df.shape[0]):
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
            if debug:
                print("Defining process: ", proc, " with id: ", cid, " and order: ", order, " in recola2.py")
            define_process_rcl(cid, proc , order)

            unselect_all_gs_powers_BornAmpl_rcl(cid)
            select_gs_power_BornAmpl_rcl(cid, 2)

            unselect_all_gs_powers_LoopAmpl_rcl(cid)
            select_gs_power_LoopAmpl_rcl(cid, 4)

    print("Generating processes in recola2.py")

    # generate it
    generate_processes_rcl()



    in_df['A0'] = 0.0
    in_df['A1'] = 0.0
    in_df['A2'] = 0.0
    for i,row in tqdm.tqdm(in_df.iterrows(),desc="Computing processes in recola2.py",total=in_df.shape[0]):
        in_ids = row[in_cols]
        out_ids = row[out_cols]
        cid = global_diags.index((*in_ids, *out_ids))
        pn_cols = sorted([col for col in in_df.columns if str(col).startswith('p_')] )
        # get the momenta of the particles
        psp = row[pn_cols].to_numpy().astype(float).reshape((int(len(pn_cols)/4),4)).tolist()

        if debug:
            print("Computing process: ", cid, " with momenta: ", psp, " and order: ", order, " in recola2.py")

        set_mu_uv_rcl(row['mu_r'])
        set_mu_ir_rcl(row['mu_r'])
        set_mu_ms_rcl(row['mu_r'])
        A0, A1, A2 = compute_process_rcl(cid, psp, order)

        alphas = get_alphas_rcl()
        alpha = get_alpha_rcl()


        if debug:
            print("alphas: ", alphas, " alpha: ", alpha)

        in_df.loc[i,('Born',)] = A0
        in_df.loc[i,('Virt',)] = A1
        in_df.loc[i,('A0',)] = A0
        in_df.loc[i,('A1',)] = A1
        in_df.loc[i,('A2',)] = A2
        in_df.loc[i,('alphas',)] = alphas
        in_df.loc[i,('alpha',)] = alpha
    return in_df


