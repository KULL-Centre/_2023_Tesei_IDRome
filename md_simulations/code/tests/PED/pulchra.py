import numpy as np
import pandas as pd
import mdtraj as md
import time
import ray
import os
import subprocess
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--num_cpus',dest='num_cpus',type=int)
parser.add_argument('--name',dest='name',type=str)
parser.add_argument('--pulchra',dest='pulchra_path',type=str)
args = parser.parse_args()

ray.init(num_cpus=args.num_cpus)

def fix_topology(prot_name):
    """
    Changes atom names to CA
    """
    trajs = []
    for replica in range(5):
        t = md.load_xtc(prot_name+f'/{replica:d}/traj.xtc',prot_name+f'/{replica:d}/top.pdb')[10:]
        cgtop = md.Topology()
        cgchain = cgtop.add_chain()
        for atom in t.top.atoms:
            cgres = cgtop.add_residue(atom.name, cgchain)
            cgtop.add_atom('CA', element=md.element.carbon, residue=cgres)
        trajs.append(md.Trajectory(t.xyz, cgtop, t.time, t.unitcell_lengths, t.unitcell_angles))
    traj = md.join(trajs)
    traj = traj.superpose(traj, frame=0)
    return traj

@ray.remote(num_cpus=1)
def run_pulchra(prot_name,pulchra_path,i,frame):
    """
    This function runs Pulchra on a single frame.
    """
    name = f'all_atom/{prot_name:s}_{i:d}.pdb'
    frame.save(name)
    FNULL = open(os.devnull, 'w')
    subprocess.run([pulchra_path,name],stdout=FNULL,stderr=FNULL)
    outname = f'all_atom/{prot_name:s}_{i:d}.rebuilt.pdb'
    trajtemp = md.load(outname)
    os.remove(name)
    os.remove(outname)
    return trajtemp.xyz

def reconstruct_pulchra(pulchra_path,prot_name,num_cpus):
    """
    This function reconstructs an all-atom trajectory from a Calpha trajectory.
    Input: trajectory nvt.xtc and topology nvt.gro file.
    n_procs: number of processors to use in parallel.
    Return: A reconstructed mdtraj trajectory.
    """
    t = fix_topology(prot_name)
    name = f'all_atom/{prot_name:s}_0.pdb'
    t[0].save(name)
    subprocess.run([pulchra_path,name])
    s = md.load_pdb(f'all_atom/{prot_name:s}_0.rebuilt.pdb')
    # n_blocks = t.n_frames // num_cpus
    xyz = np.empty((0,s.n_atoms,3))
    xyz = np.append( xyz, s.xyz )
    num_cpus = num_cpus - 1
    for j in range(1, t.n_frames, num_cpus):
        n = j+num_cpus if j+num_cpus<t.n_frames else t.n_frames
        xyz = np.append( xyz, np.vstack(ray.get([run_pulchra.remote(prot_name,pulchra_path,i,t[i]) for i in range(j,n)])) )
    allatom0 = md.Trajectory(xyz.reshape(t.n_frames,s.n_atoms,3), s.top, t.time, t.unitcell_lengths, t.unitcell_angles)
    top = md.Topology()
    chain = top.add_chain()
    for residue in allatom0.top.residues:
        res = top.add_residue(residue.name, chain, resSeq=residue.index+1)
        for atom in residue.atoms:
            top.add_atom(atom.name, element=atom.element, residue=res)
    allatom1 = md.Trajectory(allatom0.xyz, top, t.time, t.unitcell_lengths, t.unitcell_angles)
    allatom1.save_xtc(f'all_atom/{prot_name:s}.xtc')
    allatom1[0].save_pdb(f'all_atom/{prot_name:s}.pdb')
    print(prot_name,'has',allatom1.n_frames,'frames')

t0 = time.time()
reconstruct_pulchra(args.pulchra_path,args.name,args.num_cpus)
print('Timing {:.3f}'.format(time.time()-t0))
