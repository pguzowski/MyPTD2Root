#!/bin/env python

try:
    import ROOT
except ImportError, e:
    import sys
    import os
    cadbrew = None
    for p in os.environ['PATH'].split(':'):
        pp = p.split('/')
        if len(pp[-1]) == 0:
            pp.pop()
        if pp[-2] == 'CadfaelBrew' and pp[-1] == 'bin':
            cadbrew = '/'.join(pp[:-1])
            break
    if cadbrew is None:
        raise ImportError, e
    sys.path.append(os.path.join(cadbrew, 'lib','root'))
    import ROOT

import array
import math
import argparse

c = 299.792458 # speed of light in mm / ns;
m_e = 0.510998910 # electron mass in MeV

def calc_tof(ene, dEne, tlen):
    gaminv = m_e / (m_e + ene)
    beta = math.sqrt(1. - gaminv*gaminv)
    dBeta = dEne * gaminv * gaminv * gaminv / m_e / beta
    tof = tlen / beta / c
    dTof = dEne * tof * gaminv * gaminv * gaminv / beta / beta / m_e
    return (tof, dTof, beta, dBeta)


def convert(inf, outf):
    f = ROOT.TFile(inf)
    t = f.Get("PTD")
    v1 = ROOT.vector('int')()
    v1x = ROOT.vector('double')()
    v1y = ROOT.vector('double')()
    v1z = ROOT.vector('double')()
    v2 = ROOT.vector('int')()
    v2x = ROOT.vector('double')()
    v2y = ROOT.vector('double')()
    v2z = ROOT.vector('double')()
    tlen = ROOT.vector('double')()
    ene = ROOT.vector('double')()
    dEne = ROOT.vector('double')()
    time = ROOT.vector('double')()
    dTime = ROOT.vector('double')()
    Q = ROOT.vector('int')()
    np = array.array('i',[0])
    ng = array.array('i',[0])
    t.SetBranchAddress('particle.nofparticles',np)
    t.SetBranchAddress('particle.nofgammaonly',ng)
    entries = []
    for e in xrange(t.GetEntries()):
        t.GetEntry(e)
        # check number of tracks is 2 and no other gammas (unassociated calorimeter hits)
        if np[0] != 2 or ng[0] != 0:
            continue
        entries.append(e)
    
    t.SetBranchAddress('particle.vertex1_type',v1)
    t.SetBranchAddress('particle.vertex2_type',v2)
    t.SetBranchAddress('particle.vertex1_x',v1x)
    t.SetBranchAddress('particle.vertex1_y',v1y)
    t.SetBranchAddress('particle.vertex1_z',v1z)
    t.SetBranchAddress('particle.vertex2_x',v2x)
    t.SetBranchAddress('particle.vertex2_y',v2y)
    t.SetBranchAddress('particle.vertex2_z',v2z)
    t.SetBranchAddress('particle.traj_length',tlen)
    t.SetBranchAddress('particle.calo1_energy',ene)
    t.SetBranchAddress('particle.calo1_sigma_energy',dEne)
    t.SetBranchAddress('particle.calo1_time',time)
    t.SetBranchAddress('particle.calo1_sigma_time',dTime)
    t.SetBranchAddress('particle.charge',Q)

    fout = ROOT.TFile(outf,'recreate')
    tout = ROOT.TTree('t','t')

    fout_ev = array.array('i',[0]) # event
    fout_e1 = array.array('f',[0.]) # energy 1
    fout_e2 = array.array('f',[0.]) # energy 2
    fout_esum = array.array('f',[0.]) # energy sum
    fout_dv = array.array('f',[0.]) # delta vertex
    fout_pi = array.array('f',[0.]) # internal probability
    fout_pe1 = array.array('f',[0.]) # external probability 1
    fout_pe2 = array.array('f',[0.]) # external probability 2
    fout_l1 = array.array('f',[0.]) # track length 1
    fout_l2 = array.array('f',[0.]) # track length 2
    fout_lsum = array.array('f',[0.]) # track length sum
    fout_dt = array.array('f',[0.]) # time difference
    fout_q1 = array.array('i',[0]) # charge 1
    fout_q2 = array.array('i',[0]) # charge 2
    fout_qsum = array.array('i',[0]) # charge sum
    
    tout.Branch('ev',fout_ev,'ev/I')
    tout.Branch('e1',fout_e1,'e1/F')
    tout.Branch('e2',fout_e2,'e2/F')
    tout.Branch('esum',fout_esum,'esum/F')
    tout.Branch('dv',fout_dv,'dv/F')
    tout.Branch('pi',fout_pi,'pi/F')
    tout.Branch('pe1',fout_pe1,'pe1/F')
    tout.Branch('pe2',fout_pe2,'pe2/F')
    tout.Branch('l1',fout_l1,'l1/F')
    tout.Branch('l2',fout_l2,'l2/F')
    tout.Branch('lsum',fout_lsum,'lsum/F')
    tout.Branch('dt',fout_dt,'dt/F')
    tout.Branch('q1',fout_q1,'q1/I')
    tout.Branch('q2',fout_q2,'q2/I')
    tout.Branch('qsum',fout_qsum,'qsum/I')
    
    for e in entries:
        t.GetEntry(e)
        #print v1[0], v2[0], v1[1], v2[1], ene[0], ene[1]
        #each track starts on foil and ends on calo
        if not (((v1[0] == 0 and v2[0] > 1000) or (v2[0] == 0 and v1[0] > 1000)) and ((v1[1] == 0 and v2[1] > 1000) or (v2[1] == 0 and v1[1] > 1000))):
            continue

        #each track associated to calo
        if ene[0] < 0 or ene[1] < 0:
            continue

        #not same calorimeter
        if not (ene[0] > ene[1] or ene[0] < ene[1]):
            continue

        #calculate vertex difference
        if v1[0] == 0:
            vv1x = v1x[0]
            vv1y = v1y[0]
            vv1z = v1z[0]
        elif v2[0] == 0:
            vv1x = v2x[0]
            vv1y = v2y[0]
            vv1z = v2z[0]
        if v1[1] == 0:
            vv2x = v1x[1]
            vv2y = v1y[1]
            vv2z = v1z[1]
        elif v2[1] == 0:
            vv2x = v2x[1]
            vv2y = v2y[1]
            vv2z = v2z[1]
        delta_v = math.sqrt((vv1x-vv2x)*(vv1x-vv2x) + (vv1y-vv2y)*(vv1y-vv2y) + (vv1z-vv2z)*(vv1z-vv2z))
        
        #calculate internal and external probabilities
        delta_t = time[0] - time[1]

        tof_1, dTof_1, beta_1, dBeta_1 = calc_tof(ene[0], dEne[0], tlen[0])
        tof_2, dTof_2, beta_2, dBeta_2 = calc_tof(ene[1], dEne[1], tlen[1])

        delta_t_hyp = tof_1 - tof_2

        int_numer = (delta_t - delta_t_hyp)*(delta_t - delta_t_hyp)
        int_denom = dTime[0] * dTime[0] + dTime[1] * dTime[1] + dTof_1 * dTof_1 + dTof_2 * dTof_2

        int_chi2 = int_numer / int_denom
        int_prob = ROOT.TMath.Prob(int_chi2, 1)

        tof_ext_0_to_1 = (tlen[0] + tlen[1]) / (beta_2 * c)
        dTof_ext_0_to_1 = tof_ext_0_to_1 / beta_2 * dBeta_2

        ext_numer_0_to_1 = (time[1] - time[0] - tof_ext_0_to_1) * (time[1] - time[0] - tof_ext_0_to_1)
        ext_denom_0_to_1 = dTime[1] * dTime[1] + dTime[0] * dTime[0] + dTof_ext_0_to_1 * dTof_ext_0_to_1

        ext_chi2_0_to_1 = ext_numer_0_to_1 / ext_denom_0_to_1
        ext_prob_0_to_1 = ROOT.TMath.Prob(ext_chi2_0_to_1, 1)

        tof_ext_1_to_0 = (tlen[0] + tlen[1]) / (beta_1 * c)
        dTof_ext_1_to_0 = tof_ext_1_to_0 / beta_1 * dBeta_1

        ext_numer_1_to_0 = (time[0] - time[1] - tof_ext_1_to_0) * (time[0] - time[1] - tof_ext_1_to_0)
        ext_denom_1_to_0 = dTime[1] * dTime[1] + dTime[0] * dTime[0] + dTof_ext_1_to_0 * dTof_ext_1_to_0

        ext_chi2_1_to_0 = ext_numer_1_to_0 / ext_denom_1_to_0
        ext_prob_1_to_0 = ROOT.TMath.Prob(ext_chi2_1_to_0, 1)


        fout_ev[0] = e
        fout_e1[0] = ene[0]
        fout_e2[0] = ene[1]
        fout_esum[0] = ene[0]+ene[1]
        fout_dv[0] = delta_v
        fout_pi[0] = int_prob
        fout_pe1[0] = ext_prob_0_to_1
        fout_pe2[0] = ext_prob_1_to_0
        fout_l1[0] = tlen[0]
        fout_l2[0] = tlen[1]
        fout_lsum[0] = tlen[0]+tlen[1]
        fout_dt[0] = time[0] - time[1]
        fout_q1[0] = fout_q2[0] = 0
        if Q[0] == 4:
            fout_q1[0] = 1
        if Q[0] == 8:
            fout_q1[0] = -1
        if Q[1] == 4:
            fout_q2[0] = 1
        if Q[1] == 8:
            fout_q2[0] = -1
        fout_qsum[0] = fout_q1[0] + fout_q2[0]

        tout.Fill()

    tout.Write()
    fout.Close()

    
def main():
    parser = argparse.ArgumentParser(description='conversion')
    parser.add_argument('-i','--input',required=True,help='Input directory name')
    parser.add_argument('-o','--output',required=True,help='Check if this run can begin')
    args = parser.parse_args()
    convert(args.input, args.output)

if __name__ == "__main__":
    main()
