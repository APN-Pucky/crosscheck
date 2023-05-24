import subprocess
import os
from particle import Particle
import warnings
import glob
import re
import tqdm
from smpl_io import io as sio
from smpl_io import sed

ol_standalone="""
! Minimal example how to use the native interface of OpenLoops.
! It calculates the Tree and loop matrix element of the process
! u dbar -> W+ g for a random phase-space point.

! Takes a file unit, and returns the number of lines in the file.
! N.B. the file must be at the start of the file.
function count_lines(file_unit) result(output)
  integer, intent(in) :: file_unit
  integer :: output

  integer :: iostat

  output = 0
  iostat = 0
  do while (iostat==0)
    read(file_unit, *, iostat=iostat)
    if (iostat==0) then
      output = output+1
    endif
  enddo
  rewind(file_unit)
end function

program main
  use openloops
  implicit none
  integer :: id
  real(8) :: psp(0:3,5), m2l0, m2l1(0:2), acc, sqrt_s = 1000
  integer              :: in_unit,lines,i,count_lines
  character(100) :: key
  real(8) :: val

  !call set_parameter("order_ew", 1) ! coupling order
  !call set_parameter("ew_scheme", 2) ! coupling order
  open (newunit=in_unit, file = "params.input", action = "read", status = "old")
  lines = count_lines(in_unit)
  do i=1,lines
        read(in_unit,*) key, val
        call set_parameter(key,val)
  enddo
  id = register_process("REPLACETHISBYTHEPROCESS", 11) ! register loop process d dbar --> Z u ubar
  call start() ! initialise

  open (newunit=in_unit, file = "PS.input", action = "read", status = "old")
  lines = count_lines(in_unit)
  !print *,lines
  do i=1,lines
        read(in_unit,*) psp(:,i)
  enddo

  print *,psp(:,1)
  !call exit(0)

  call evaluate_loop(id, psp, m2l0, m2l1, acc) ! calculate tree and loop matrix elements
  call finish()

  print *, "  u ubar --> A c cbar for a RAMBO phase space point"
  print *, "  tree                     loop ep^0                  loop ep^-1               loop ep^-2                 accuracy"
  print *, "born = ", m2l0
  print *, "finite = ", m2l1(0)
  print *, "1eps = ", m2l1(1)
  print *, "2eps = ", m2l1(2)
  print *, "acc = ", acc
end program main
"""

ol_params="""order_ew 1
ew_scheme 2
mass(23) 91.2
mass(1) 0.0d0
mass(2) 0.0d0
mass(3) 0.0d0
mass(4) 0.0d0
mass(5) 0.0d0
"""

def compute(in_df, debug=False,n_jobs=-1):
    global_diags = []
    procs=[]
    in_cols = sorted([col for col in in_df.columns if str(col).startswith('in_')] )
    out_cols = sorted([col for col in in_df.columns if str(col).startswith('out_')])

    for i,row in tqdm.tqdm(in_df.iterrows(),desc="Generating processes in openloops.py",total=in_df.shape[0]):
        # get the ids of the particles
        in_ids = [row[col] for col in in_cols]
        out_ids = [row[col] for col in out_cols]
        # do not append to global_diags if it is already there
        if not (*in_ids, *out_ids) in global_diags:
            global_diags.append((*in_ids, *out_ids))
            cid = len(global_diags) - 1
            if not os.path.exists(f"/tmp/ol_{cid}/"):
                os.mkdir(f"/tmp/ol_{cid}/")
            #gfortran -o OL_fortran -Wl,-rpath=../lib -I../lib_src/openloops/mod -L../lib -lopenloops OL_fortran.f90
            with open(f"/tmp/ol_{cid}/OL_fortran.f90",'w') as f:
                f.write(ol_standalone.replace("REPLACETHISBYTHEPROCESS",f"{' '.join([str(i) for i in in_ids])} -> {' '.join([str(i) for i in out_ids])}") )
            if not n_jobs < 0 and len(procs) > n_jobs:
                procs[len(procs)-n_jobs].wait()
            cmp = ["gfortran", "-o", "OL_fortran", "-Wl,-rpath=/opt/OpenLoops-2.1.2/proclib","-I/usr/include/","-lqcdloop","-lcuttools","-lopenloops","OL_fortran.f90"]
            if debug:
                print(" ".join(cmp))
            procs.append(subprocess.Popen(cmp,
                                          stdout=subprocess.DEVNULL if not debug else None, 
                                          cwd=f"/tmp/ol_{cid}/"))
        
    for proc in procs:
        proc.wait()


    for i,row in tqdm.tqdm(in_df.iterrows(),desc="Computing processes in openloops.py",total=in_df.shape[0]):
        in_ids = row[in_cols]
        out_ids = row[out_cols]
        cid = global_diags.index((*in_ids, *out_ids))
        sio.write(f"/tmp/ol_{cid}/params.input",
                  ol_params + 
                  f"mu {row['mu_r']}\n"+
                  f"alpha_s {row['alphas']}\n"+
                  f"alpha_qed {row['alpha']}\n"
                  )
        pn_cols = sorted([col for col in in_df.columns if str(col).startswith('p_')] )
        psp = row[pn_cols].to_numpy().astype(float).reshape((int(len(pn_cols)/4),4)).tolist()
        with open(f"/tmp/ol_{cid}/PS.input",'w') as f:
            for ps in psp:
                print(*ps, file=f)
        output = subprocess.check_output([f"/tmp/ol_{cid}/OL_fortran"], cwd=f"/tmp/ol_{cid}/").decode('utf-8')
        if debug:
            print(output)
        born = re.search(r"born\s*=\s*(.*)", output).group(1)
        finite = re.search(r"finite\s*=\s*(.*)", output).group(1)
        single = re.search(r"1eps\s*=\s*(.*)", output).group(1)
        double = re.search(r"2eps\s*=\s*(.*)", output).group(1)
        in_df.loc[i,('Born',)] = float(born)
        in_df.loc[i,('Virt',)] = float(finite)
        in_df.loc[i,('Finite',)] = float(finite)
        in_df.loc[i,('Single Pole',)] = float(single)
        in_df.loc[i,('Double Pole',)] = float(double)

                      