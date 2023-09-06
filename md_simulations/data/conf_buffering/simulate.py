import hoomd
import hoomd.md
import time
import os
import sys
import itertools
import pandas as pd
import numpy as np
import mdtraj as md
from hoomd import azplugins
from argparse import ArgumentParser
from scipy.optimize import curve_fit
from itertools import combinations
from numpy import linalg
from sklearn.covariance import LedoitWolf
from scipy import constants
sys.path.append('BLOCKING')
from main import BlockAnalysis

parser = ArgumentParser()
parser.add_argument('--seq_name',nargs='?',required=True)
parser.add_argument('--path',nargs='?',required=True)
args = parser.parse_args()

print(hoomd.__file__)

def xy_spiral_array(n, delta=0, arc=.38, separation=.7):
    """
    create points on an Archimedes' spiral
    with `arc` giving the length of arc between two points
    and `separation` giving the distance between consecutive
    turnings
    """
    def p2c(r, phi):
        """
        polar to cartesian
        """
        return (r * np.cos(phi), r * np.sin(phi))
    r = arc
    b = separation / (2 * np.pi)
    phi = float(r) / b
    coords = []
    for i in range(n):
        coords.append(list(p2c(r, phi))+[0])
        phi += float(arc) / r
        r = b * phi
    return np.array(coords)+delta

def genTop(residues,fasta,path,L):
    N_res = len(fasta)
    top = md.Topology()
    chain = top.add_chain()
    for resname in fasta:
        residue = top.add_residue(residues.loc[resname,'three'], chain)
        top.add_atom(residues.loc[resname,'three'], element=md.element.carbon, residue=residue)
    for i in range(N_res-1):
        top.add_bond(top.atom(i),top.atom(i+1))
    pos = [[0,0,(i-N_res/2.)*.38] for i in range(N_res)]
    t = md.Trajectory(np.array(pos).reshape(N_res,3), top, 0, [L,L,L], [90,90,90])
    t.save_pdb(path+'/top.pdb')

def genParams(r,seq,temp,ionic):
    RT = 8.3145*temp*1e-3
    pH = 7.4
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(temp)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/RT
    # Calculate the inverse of the Debye length
    yukawa_kappa = np.sqrt(8*np.pi*lB*ionic*6.022/10)
    fasta = list(seq.fasta)
    #r.loc['H','q'] = 1. / ( 1 + 10**(pH-6) )
    if seq['first'] == 1:
        r.loc['X'] = r.loc[fasta[0]]
        r.loc['X','q'] = r.loc[fasta[0],'q'] + 1.
        r.loc['X','MW'] = r.loc[fasta[0],'MW'] + 2.
        fasta[0] = 'X'
    if seq['last'] == seq.N_FL:
        r.loc['Z'] = r.loc[fasta[-1]]
        r.loc['Z','q'] = r.loc[fasta[-1],'q'] - 1.
        r.loc['Z','MW'] = r.loc[fasta[-1],'MW'] + 16.
        fasta[-1] = 'Z'
    # Calculate the prefactor for the Yukawa potential
    qq = pd.DataFrame(r.q.values*r.q.values.reshape(-1,1),index=r.q.index,columns=r.q.index)
    yukawa_eps = qq*lB*RT
    types = list(np.unique(fasta))
    pairs = np.array(list(itertools.combinations_with_replacement(types,2)))
    return yukawa_kappa, yukawa_eps, types, pairs, fasta, r

#energy functions
HALR = lambda r,s,l : 4*0.8368*l*((s/r)**12-(s/r)**6)
HASR = lambda r,s,l : 4*0.8368*((s/r)**12-(s/r)**6)+0.8368*(1-l)
HA = lambda r,s,l : np.where(r<2**(1/6)*s, HASR(r,s,l), HALR(r,s,l))
HASP = lambda r,s,l,rc : np.where(r<rc, HA(r,s,l)-HA(rc,s,l), 0)

#force functions
LJ_F = lambda r,s,rvec : -6*4*0.8368*(2*s**12/r**14-s**6/r**8)*rvec
HA_F = lambda r,s,l,rvec : np.where(r<2**(1/6)*s, LJ_F(r,s,rvec), l*LJ_F(r,s,rvec))
HASP_F = lambda r,s,l,rvec,rc : np.where(r<rc, HA_F(r,s,l,rvec), 0)

DH_F = lambda r,yukawa_eps,yukawa_kappa,rvec : -yukawa_eps*np.exp(-r*yukawa_kappa)*(1+r*yukawa_kappa)/r**3*rvec
DHSP_F = lambda r,yukawa_eps,yukawa_kappa,rvec,rc : np.where(r<rc, DH_F(r,yukawa_eps,yukawa_kappa,rvec), 0)

def calcEnergyMap(t,df,prot,rc):
    indices = t.top.select_pairs('all','all')
    mask = np.abs(indices[:,0]-indices[:,1])>1 #exclude >1, was used to exclude bonded pairs
    indices = indices[mask]
    dvec = md.compute_displacements(t,indices) #vector distances between pairs for each frame
    d = linalg.norm(dvec,axis=2)
    pairs = np.array(list(combinations(list(prot.fasta),2)))
    pairs = pairs[mask]
    sigmas = 0.5*(df.loc[pairs[:,0]].sigmas.values+df.loc[pairs[:,1]].sigmas.values)
    lambdas = 0.5*(df.loc[pairs[:,0]].lambdas.values+df.loc[pairs[:,1]].lambdas.values)
    temp = 310
    ionic = 0.15
    RT = 8.3145*temp*1e-3
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(temp)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/RT
    # Calculate the inverse of the Debye length
    yukawa_kappa = np.sqrt(8*np.pi*lB*ionic*6.022/10)
    qq = df.loc[pairs[:,0]].q.values*df.loc[pairs[:,1]].q.values
    yukawa_eps = qq*lB*RT
    emap = np.zeros(pairs.shape[0])
    forces = np.zeros((t.n_frames,t.n_atoms,3))
    dstd = np.std(d,axis=0)
    for j,(r,rvec) in enumerate(zip(np.split(d,20,axis=0),np.split(dvec,20,axis=0))):
        emap += np.nansum(HASP(r,sigmas[np.newaxis,:],lambdas[np.newaxis,:],rc),axis=0)
        for i in range(t.n_atoms):
            ndx_pairs = np.any(indices==i,axis=1)
            forces[j*r.shape[0]:(j+1)*r.shape[0],i] = np.nansum(HASP_F(r[:,ndx_pairs,np.newaxis],sigmas[np.newaxis,ndx_pairs,np.newaxis],
                            lambdas[np.newaxis,ndx_pairs,np.newaxis],rvec[:,ndx_pairs],rc),axis=1)
            forces[j*r.shape[0]:(j+1)*r.shape[0],i] += np.nansum(DHSP_F(r[:,ndx_pairs,np.newaxis],
                            yukawa_eps[np.newaxis,ndx_pairs,np.newaxis],yukawa_kappa,rvec[:,ndx_pairs],4),axis=1)
    return indices, emap/d.shape[0], forces, dstd

def calcRg(t,residues,seq):
    fasta = list(seq.fasta)
    masses = residues.loc[fasta,'MW'].values
    # calculate the center of mass
    cm = np.sum(t.xyz*masses[np.newaxis,:,np.newaxis],axis=1)/masses.sum()
    # calculate residue-cm distances
    si = np.linalg.norm(t.xyz - cm[:,np.newaxis,:],axis=2)
    # calculate rg
    rgarray = np.sqrt(np.sum(si**2*masses,axis=1)/masses.sum())
    return rgarray

def error_ratio(v1,v2,e1,e2):
    ratio = v1/v2
    return ratio*np.sqrt((e1/v1)**2+(e2/v2)**2)

def calcRgTensor(t,residues,seq,forces):
    fasta = list(seq.fasta)
    masses = residues.loc[fasta,'MW'].values[np.newaxis,:,np.newaxis]

    prefactor = constants.Boltzmann*310/constants.hbar**2
    forces = forces/np.sqrt(masses/1e3/constants.Avogadro)
    forces = forces.reshape(t.n_frames,-1)*1e3/constants.Avogadro*1e9
    sigma = LedoitWolf().fit(forces).covariance_
    eigenvalues = linalg.eigvalsh(sigma)
    kT_over_hbar_omega = constants.Boltzmann*310*np.sqrt(prefactor/eigenvalues)
    SPR = np.sum(np.log(kT_over_hbar_omega+1))/seq.N # R

    # calculate the center of mass
    cm = np.sum(t.xyz*masses,axis=1)/masses.sum()
    # calculate residue-cm distances
    si = t.xyz - cm[:,np.newaxis,:]
    q = np.einsum('jim,jin->jmn', si*masses,si)/masses.sum()
    trace_q = np.trace(q,axis1=1,axis2=2)
    # calculate rg
    rgarray = np.sqrt(trace_q)
    # calculate traceless matrix
    mean_trace = np.trace(q,axis1=1,axis2=2)/3
    q_hat = q - mean_trace.reshape(-1,1,1)*np.identity(3).reshape(-1,3,3)
    # calculate asphericity
    Delta_array = 3/2*np.trace(q_hat**2,axis1=1,axis2=2)/(trace_q**2)
    # calculate oblateness
    S_array = 27*linalg.det(q_hat)/(trace_q**3)
    # calculate ensemble averages
    block_tr_q_hat_2 = BlockAnalysis(np.trace(q_hat**2,axis1=1,axis2=2), multi=1)
    block_tr_q_hat_2.SEM()
    block_tr_q_2 = BlockAnalysis(trace_q**2, multi=1)
    block_tr_q_2.SEM()
    block_det_q_hat = BlockAnalysis(linalg.det(q_hat), multi=1)
    block_det_q_hat.SEM()
    block_tr_q_3 = BlockAnalysis(trace_q**3, multi=1)
    block_tr_q_3.SEM()
    Delta = 3/2*block_tr_q_hat_2.av/block_tr_q_2.av
    S = 27*block_det_q_hat.av/block_tr_q_3.av
    Delta_err = 3/2*error_ratio(block_tr_q_hat_2.av,block_tr_q_2.av,block_tr_q_hat_2.sem,block_tr_q_2.sem)
    S_err = 27*error_ratio(block_det_q_hat.av,block_tr_q_3.av,block_det_q_hat.sem,block_tr_q_3.sem)
    return rgarray, Delta_array, S_array, Delta, S, Delta_err, S_err, SPR

def calcRs(traj):
    pairs = traj.top.select_pairs('all','all')
    d = md.compute_distances(traj,pairs)
    nres = traj.n_atoms
    ij = np.arange(2,nres,1)
    diff = [x[1]-x[0] for x in pairs]
    dij = np.empty(0)
    for i in ij:
        dij = np.append(dij,np.sqrt((d[:,diff==i]**2).mean().mean()))
    return ij,dij,np.mean(1/d,axis=1)

def analyse(residues,path,seq):
    top = md.Topology()
    chain = top.add_chain()
    if os.path.exists(path+'/traj.gsd'):
        traj = md.load(path+'/traj.gsd')
    else:
        traj = md.load_xtc(path+'/traj.xtc',top=path+'/top.pdb')
    N_res = traj.n_atoms
    fasta = list(seq.fasta)
    #fixing the trajectory to the middle of the box
    for resname in fasta:
        residue = top.add_residue(residues.loc[resname,'three'],chain)
        top.add_atom(residues.loc[resname,'three'], element=md.element.carbon, residue=residue)
    for i in range(N_res-1):
        top.add_bond(top.atom(i),top.atom(i+1))
    traj.top = top
    traj = traj.image_molecules(inplace=False, anchor_molecules=[set(traj.top.atoms)],make_whole=True)
    print('Number of frames: {:d}'.format(traj.n_frames))
    if os.path.exists(path+'/traj.gsd'):
        traj[-1].save_pdb(path+'/top.pdb')
        traj.save_xtc(path+'/traj.xtc')
        os.remove(path+'/traj.gsd')
    #skip first 10 frames
    traj = traj[10:]
    #energy maps
    df_emap = pd.DataFrame(index=range(traj.n_atoms),columns=range(traj.n_atoms),dtype=float)
    df_dstd = pd.DataFrame(index=range(traj.n_atoms),columns=range(traj.n_atoms),dtype=float)
    pairs, emap, forces, dstd = calcEnergyMap(traj,residues,seq,2.0)
    for k,(i,j) in enumerate(pairs):
        df_emap.loc[i,j] = emap[k]
        df_emap.loc[j,i] = emap[k]
        df_dstd.loc[i,j] = dstd[k]
        df_dstd.loc[j,i] = dstd[k]
    df_analysis = pd.DataFrame(index=['Rg','ete','rh','nu','R0','ete2_Rg2','Delta','S','SPR'],columns=['value','error'])
    #rg
    rgarray, Delta_array, S_array, Delta, S, Delta_err, S_err, SPR = calcRgTensor(traj,residues,seq,forces)
    df_analysis.loc['Delta','value'] = Delta
    df_analysis.loc['Delta','error'] = Delta_err
    df_analysis.loc['S','value'] = S
    df_analysis.loc['S','error'] = S_err
    df_analysis.loc['SPR','value'] = SPR
    df_analysis.loc['SPR','error'] = 0
    np.save(path+'/rg.npy',rgarray)
    np.save(path+'/Delta.npy',Delta_array)
    np.save(path+'/S.npy',S_array)
    block_rg = BlockAnalysis(rgarray, multi=1)
    block_rg.SEM()
    df_analysis.loc['Rg','value'] = block_rg.av
    df_analysis.loc['Rg','error'] = block_rg.sem
    #ete
    ete = md.compute_distances(traj,atom_pairs=[[0,N_res-1]]).flatten()
    np.save(path+'/ete.npy',ete)
    block_ete = BlockAnalysis(ete, multi=1)
    block_ete.SEM()
    df_analysis.loc['ete','value'] = block_ete.av
    df_analysis.loc['ete','error'] = block_ete.sem
    block_ete2 = BlockAnalysis(np.power(ete,2), multi=1)
    block_ete2.SEM()
    block_rg2 = BlockAnalysis(np.power(rgarray,2), multi=1)
    block_rg2.SEM()
    ete2 = block_ete2.av
    rg2 = block_rg2.av
    ete2_e = block_ete2.sem
    rg2_e = block_rg2.sem
    df_analysis.loc['ete2_Rg2','value'] = ete2 / rg2
    df_analysis.loc['ete2_Rg2','error'] = error_ratio(ete2,rg2,ete2_e,rg2_e)
    #nonlinear scaling exponent
    f = lambda x,R0,v : R0*np.power(x,v)
    ij,dij,invrij = calcRs(traj)
    block_invrij = BlockAnalysis(invrij, multi=1)
    block_invrij.SEM()
    df_analysis.loc['rh','value'] = 1/(1-1/N_res)/block_invrij.av
    df_analysis.loc['rh','error'] = block_invrij.sem/(1-1/N_res)/block_invrij.av/block_invrij.av
    np.save(path+'/rs.npy',dij)
    popt, pcov = curve_fit(f,ij[ij>5],dij[ij>5],p0=[.4,.5])
    df_analysis.loc['nu','value'] = popt[1]
    df_analysis.loc['nu','error'] = pcov[1,1]**0.5
    df_analysis.loc['R0','value'] = popt[0]
    df_analysis.loc['R0','error'] = pcov[0,0]**0.5
    return df_emap,df_dstd,df_analysis

def simulate(residues,sequences,seq_name,path):
    seq = sequences.loc[seq_name]
    N_res = seq.N

    if N_res > 500:
        hoomd.context.initialize("--mode=gpu");
        L = 200
    else:
        hoomd.context.initialize("--mode=cpu");
        L = (N_res-1)*0.38+4

    hoomd.option.set_notice_level(1)
    hoomd.util.quiet_status()

    lj_eps = 4.184*.2
    temp = 310
    ionic_strength = 0.15 # M
    RT = 8.3145*temp*1e-3

    yukawa_kappa, yukawa_eps, types, pairs, fasta, residues = genParams(residues,seq,temp,ionic_strength)

    sigmamap = pd.DataFrame((residues.sigmas.values+residues.sigmas.values.reshape(-1,1))/2,
                            index=residues.sigmas.index,columns=residues.sigmas.index)
    lambdamap = pd.DataFrame((residues.lambdas.values+residues.lambdas.values.reshape(-1,1))/2,
                            index=residues.lambdas.index,columns=residues.lambdas.index)

    N_save = 7000 if N_res < 150 else int(np.ceil(3e-4*N_res**2)*1000)
    N_steps = 1010*N_save

    genTop(residues,fasta,path,L)

    snapshot = hoomd.data.make_snapshot(N=N_res,
                                box=hoomd.data.boxdim(Lx=L, Ly=L, Lz=L),
                                particle_types=types,
                                bond_types=['polymer']);

    snapshot.bonds.resize(N_res-1);

    if N_res > 500:
        snapshot.particles.position[:] = xy_spiral_array(N_res)
    else:
        snapshot.particles.position[:] = [[0,0,(i-N_res/2.)*.38] for i in range(N_res)]
    snapshot.particles.typeid[:] = [types.index(a) for a in fasta]
    snapshot.particles.mass[:] = [residues.loc[a].MW for a in fasta]

    snapshot.bonds.group[:] = [[i,i+1] for i in range(N_res-1)];
    snapshot.bonds.typeid[:] = [0] * (N_res-1)

    hoomd.init.read_snapshot(snapshot);

    hb = hoomd.md.bond.harmonic();
    hb.bond_coeff.set('polymer', k=8033.0, r0=0.38);

    nl = hoomd.md.nlist.cell();

    ah = azplugins.pair.ashbaugh(r_cut=2.0, nlist=nl)
    yukawa = hoomd.md.pair.yukawa(r_cut=4.0, nlist=nl)
    for a,b in pairs:
        ah.pair_coeff.set(a, b, lam=lambdamap.loc[a,b], epsilon=lj_eps, sigma=sigmamap.loc[a,b], r_cut=2.0)
        yukawa.pair_coeff.set(a, b, epsilon=yukawa_eps.loc[a,b], kappa=yukawa_kappa, r_cut=4.)

    ah.set_params(mode='shift')
    yukawa.set_params(mode='shift')
    nl.reset_exclusions(exclusions = ['bond'])

    integrator_mode = hoomd.md.integrate.mode_standard(dt=0.005);
    integrator = hoomd.md.integrate.langevin(group=hoomd.group.all(),kT=RT,seed=np.random.randint(100));

    for a in types:
        integrator.set_gamma(a, residues.loc[a].MW/100)

    hoomd.run(10000)
    integrator.disable()

    integrator_mode = hoomd.md.integrate.mode_standard(dt=0.01);
    integrator = hoomd.md.integrate.langevin(group=hoomd.group.all(),kT=RT,seed=np.random.randint(100));

    for a in types:
        integrator.set_gamma(a, residues.loc[a].MW/100)

    hoomd.dump.gsd(filename=path+'/traj.gsd', period=N_save, group=hoomd.group.all(), overwrite=True);

    hoomd.run(N_steps)

    hoomd.dump.gsd(filename=path+'/restart.gsd', group=hoomd.group.all(), truncate=True, period=None, phase=0)

residues = pd.read_csv('residues.csv').set_index('one',drop=False)

sequences = pd.read_csv('conf_buffering_seq.csv',index_col=0)

t0 = time.time()
simulate(residues,sequences,args.seq_name,args.path)
df_emap,df_dstd,df_analysis = analyse(residues,args.path,sequences.loc[args.seq_name])
df_analysis.to_csv(args.path+'/analysis.csv')
df_emap.to_csv(args.path+'/emap.csv')
df_dstd.to_csv(args.path+'/stdev.csv')
print('Timing sim and analysis {:.3f}'.format((time.time()-t0)/3600))
