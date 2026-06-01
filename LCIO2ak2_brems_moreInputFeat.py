from __future__ import print_function

import os,pprint,glob,json

from pyLCIO.io import LcioReader
from pyLCIO.UTIL import LCRelationNavigator

#import ROOT
import numpy as np
#import re

import awkward as ak
#import awkward0
#import pyarrow as pa
#import pyarrow.parquet as pq

import sys
import math
#ROOT.gROOT.SetBatch()

import h5py

def printMc(mc):
    print("ID:", mc.id(), "PDG:", mc.getPDG(), "Momentum:(", f'{mc.getMomentum()[0]:.3f}',f'{mc.getMomentum()[1]:.3f}',f'{mc.getMomentum()[2]:.3f}',") Vertex:(", f'{mc.getVertex()[0]:.3f}',f'{mc.getVertex()[1]:.3f}',f'{mc.getVertex()[2]:.3f}',")")

def isDeltaRay(mc):
    if(mc.getPDG() != 22): return 0
    pmc = mc.getParents()[0]
    if(pmc.getPDG() == 111): return 0

    print("Photon found.")
    printMc(mc)

    print("Parent:")
    printMc(pmc)

    print("Daughters:")
    for dmc in pmc.getDaughters():
        printMc(dmc)

    return 0



skimmedIds = []     # skimmed Id

def isBrems(mc):
    result = True if mc.getPDG() == 22 and mc.vertexIsNotEndpointOfParent() else False
    if mc.getPDG() in [11, -11]:
        for daughter in mc.getDaughters():
            if (daughter.id() in skimmedIds and daughter.getPDG() == 22 and daughter.vertexIsNotEndpointOfParent()):
                result = True
    return result

def isConverted(mc):
    if not mc.getParents(): return False
    if((mc.getPDG() in [11, -11]) and mc.getParents()[0].getPDG()==22):
        pdg = 11 if mc.getPDG() == -11 else -11
        for daughter in mc.getParents()[0].getDaughters():
            if mc.id() == daughter.id(): continue
            if daughter.getPDG() == pdg:
                return True
    return False

def returnDictMCP(mp):
    dictMCPId = 0
    isCon = isConverted(mp)
    isBre = isBrems(mp)
    if isCon and isBre: dictMCPId = 3
    elif isCon: dictMCPId = 1
    elif isBre: dictMCPId = 2
    return dictMCPId

def checkPi02gamma(mc):
    result = False
    if mc.getPDG() == 22 and mc.getParents()[0].getPDG() == 111:
        if mc.getParents()[0].getDaughters().size() != 2:
            return False
        result = True
        for daughter in mc.getParents()[0].getDaughters():
            result = result and daughter.getPDG() == 22
    return result

def checkBrems(mc):
    result = False
    if mc.getPDG() == 22 and mc.vertexIsNotEndpointOfParent() and mc.getParents()[0].getPDG() in [11, -11]:
        result = True
    return result

def checkDistance(pos1, pos2):
    if pos1 is None or pos2 is None:
        return False
    distance = np.linalg.norm(pos1 - pos2)
    print(pos1, pos2, distance)

    return distance < 20

B_FIELD_TESLA = 3.5

def get_track_momentum_at_point(track, x, y, z):
    """
    trackの指定点に最も近いTrackStateから運動量を計算。
    p = 0.3 * B / |omega| (GeV/c), omegaは1/mm。
    """
    state = track.getClosestTrackState(x, y, z)
    if state is None:
        return None
    p = -.001 / state.getOmega() * 0.3 * B_FIELD_TESLA
    px = p * math.cos(state.getPhi())
    py = p * math.sin(state.getPhi())
    pz = p * state.getTanLambda()
    return (px, py, pz)

def find_track_for_mcp_at_brems(mc, event, TrackList, RelTrackList, dictSkimmed, colnames):
    """
    MC粒子に紐づくtrackのうち、brems後のセグメント(最外周)を返す。
    brems点に最も近いtrack stateを持つtrackを優先。
    """
    brems_daughters = [d for d in mc.getDaughters() if d.getPDG() == 22 and d.vertexIsNotEndpointOfParent()]
    if not brems_daughters:
        return None
    brems_vtx = brems_daughters[0].getVertex()
    bx, by, bz = brems_vtx[0], brems_vtx[1], brems_vtx[2]

    best_track = None
    best_rin = -1.0

    for colname in TrackList:
        if colname not in colnames: continue
        col = event.getCollection(colname)
        for track in col:
            rin = track.getRadiusOfInnermostHit()
            if rin <= 0:
                # fallback: track stateの参照点からinnermost(最小半径)を計算
                rin = 1e30
                for ts in track.getTrackStates():
                    ref = ts.getReferencePoint()
                    r = math.sqrt(ref[0]**2 + ref[1]**2 + ref[2]**2)
                    if r < rin:
                        rin = r
                if rin >= 1e30:
                    rin = -1.0
            if rin < 0:
                continue

            for relcolname in RelTrackList:
                if relcolname not in colnames: continue
                relcol = event.getCollection(relcolname)
                nav = LCRelationNavigator(relcol)
                for rel, w in zip(nav.getRelatedToObjects(track), nav.getRelatedToWeights(track)):
                    if w <= 0.5:
                        continue
                    mapped_mcp = dictSkimmed[rel.id()] if rel.id() in dictSkimmed else rel
                    if mapped_mcp.id() != mc.id():
                        continue
                    if rin > best_rin:
                        best_rin = rin
                        best_track = track
                    break

    return best_track

def _get_brems_photons(mc):
    """
    複数回bremsしてもelectronは同じ粒子のまま。brems光子はすべて直接の娘。
    """
    return [d for d in mc.getDaughters() if d.getPDG() == 22 and d.vertexIsNotEndpointOfParent()]

def _is_photon_merged_to_mc(photon, mc, dictSkimmed):
    """
    photonがこのmcにmergeされているか。
    dictSkimmed[photon_id] == mc ならmerge済み（差し引かない）。
    """
    if photon.id() not in dictSkimmed:
        return False
    return dictSkimmed[photon.id()].id() == mc.id()

def get_brems_corrected_energy_momentum(mc, event, TrackList, RelTrackList, dictSkimmed, colnames):
    """
    bremsを起こしているe-/e+について補正する。
    - Energy: mc.getEnergy()はbrems光子分も含むため、光子分を差し引く。
    - merge済みの光子(photon_n, photon_(n-1)など)は差し引かない。photon_1〜photon_(n-2)のみ差し引く。
    - Momentum: trackのbrems点での運動量を使用（磁場補正）。merge済み光子の運動量は引かない。
    """
    energy = mc.getEnergy()
    px, py, pz = mc.getMomentum()[0], mc.getMomentum()[1], mc.getMomentum()[2]
    if mc.getPDG() not in [11, -11]:
        return energy, (px, py, pz)

    all_brems = _get_brems_photons(mc)
    # merge済みの光子は差し引かない
    to_subtract = [d for d in all_brems if not _is_photon_merged_to_mc(d, mc, dictSkimmed)]
    if not to_subtract:
        return energy, (px, py, pz)

    for d in to_subtract:
        energy -= d.getEnergy()

    track = find_track_for_mcp_at_brems(mc, event, TrackList, RelTrackList, dictSkimmed, colnames)
    if track is not None and all_brems:
        brems_vtx = max(all_brems, key=lambda d: math.sqrt(d.getVertex()[0]**2 + d.getVertex()[1]**2 + d.getVertex()[2]**2)).getVertex()
        mom_at_brems = get_track_momentum_at_point(track, brems_vtx[0], brems_vtx[1], brems_vtx[2])
        if mom_at_brems is not None:
            px, py, pz = mom_at_brems[0], mom_at_brems[1], mom_at_brems[2]
        else:
            px, py, pz = mc.getMomentum()[0], mc.getMomentum()[1], mc.getMomentum()[2]
            for d in to_subtract:
                px -= d.getMomentum()[0]
                py -= d.getMomentum()[1]
                pz -= d.getMomentum()[2]
    else:
        px, py, pz = mc.getMomentum()[0], mc.getMomentum()[1], mc.getMomentum()[2]
        for d in to_subtract:
            px -= d.getMomentum()[0]
            py -= d.getMomentum()[1]
            pz -= d.getMomentum()[2]

    return energy, (px, py, pz)

def get_entry_point_from_simhits(mcp, event, calo_hits):
    """
    mcp: EVENT.MCParticle
    calo_hits: LCCollection of SimCalorimeterHit
    return: np.array([x,y,z]) or None
    """
    vtx = np.array([mcp.getVertex()[0], mcp.getVertex()[1], mcp.getVertex()[2]])
    mom = np.array([mcp.getMomentum()[0], mcp.getMomentum()[1], mcp.getMomentum()[2]])
    pnorm = np.linalg.norm(mom)
    if pnorm == 0:
        return None
    phat = mom / pnorm

    best_t = 1e30
    best_pos = None
    # for colname in calo_hits:
    #     col = event.getCollection(colname)
    #     # print("col : ", col)
    #     for simhit in col:
    #         for i in range(simhit.getNMCContributions()):
    #             if simhit.getParticleCont(i) is mcp:
    #                 pos = np.array([simhit.getPosition()[0], simhit.getPosition()[1], simhit.getPosition()[2]])
    #                 d = pos - vtx
    #                 t = np.dot(d, phat)
    #                 if t > 0 and t < best_t:
    #                     best_t = t
    #                     best_pos = pos
    for colname in calo_hits:
        col = event.getCollection(colname)
        # print("col : ", col)
        for simhit in col:
            for i in range(simhit.getNMCContributions()):
                if simhit.getParticleCont(i) is mcp:
                    pos = np.array([simhit.getPosition()[0], simhit.getPosition()[1], simhit.getPosition()[2]])
                    d_ = pos - vtx
                    t = np.dot(d_, phat)
                    d = pos - vtx + t * phat
                    if t > 0 and t < best_t:
                        best_t = t
                        best_pos = pos
    return best_pos

def get_calo_entry_from_mcp_via_track(mcp, event, TrkList, RelTrkList):
    """
    Parameters
    ----------
    mcp : EVENT.MCParticle
        対象の MCParticle
    tracks : LCCollection (Track)
        MarlinTrkTracks
    
    Returns
    -------
    np.array([x,y,z]) or None
        calorimeter 入射位置
    """
    result = None

    for colname in TrkList:
        col = event.getCollection(colname)
        for track in col:
            nstates = track.getTrackStates().size()
            if nstates > 0:
                state = track.getTrackState(4) # at calorimeter
                point = state.getReferencePoint()

                for relcolname in RelTrkList:
                    relcol = event.getCollection(relcolname)
                    nav = LCRelationNavigator(relcol)
                    for rel,w in zip(nav.getRelatedToObjects(track),nav.getRelatedToWeights(track)):
                        if w > 0.5:
                            if rel is mcp:
                                result = np.array(point[0], point[1], point[2])
    return result



#===========================================================================================
def makeAk(filename, outfilename, maxread, skip):

    reader = LcioReader.LcioReader(filename)

    #EcalParaListDict=makeEmpty(["ID","MC_PDG","MC_Status","Ecal_E","Ecal_x","Ecal_y","Ecal_z","Ecal_t","HitTag"])

    TrackList=["MarlinTrkTracks"]
    RelTrackList=["MarlinTrkTracksMCTruthLink"]
    SimHitNameList=["EcalBarrelCollection","EcalEndcapRingCollection","EcalEndcapsCollection","HcalBarrelRegCollection","HcalEndcapRingCollection","HcalEndcapsCollection"]
    CaloHitNameList=["EcalBarrelCollectionRec","EcalEndcapRingCollectionRec","EcalEndcapsCollectionRec","HcalBarrelCollectionRec","HcalEndcapRingCollectionRec","HcalEndcapsCollectionRec"]
    RelHitList=["EcalBarrelRelationsSimRec","EcalEndcapRingRelationsSimRec","EcalEndcapsRelationsSimRec","HcalBarrelRelationsSimRec","HcalEndcapRingRelationsSimRec","HcalEndcapsRelationsSimRec"]

    PandoraList=["PandoraPFOs"]
    PandoraClusterList=["PandoraClusters"]
    RelPandoraList=["RecoMCTruthLink"]

    MCParticleList=["MCParticle"]
    SkimmedList=["MCParticlesSkimmed"]

    nbevents = reader.getNumberOfEvents()
    print("# Events in the file:", nbevents, ", # events to read:", maxread)
    #if ((maxread > 0) and (nbevents > maxread)): nbevents = maxread
    nread = maxread if maxread > 0 and nbevents > maxread else nbevents
    #print(nread, nbevents, maxread)

    # =======================================================================
    # Read events in the file, and fill Ntuple
    # =======================================================================

    array_exist = 0

    abnormal_virtual_hit = 0
    notStopped = 0
    decayedInTracker = 0
    npi = 0
    nk = 0
    npi0 = 0
    n_decayedpi0 = 0
    n_decayedpi0_near = 0
    n_brems = 0
    n_brems_near = 0

    for idx,event in enumerate(reader):
        b_feat = ak.ArrayBuilder()
        b_label = ak.ArrayBuilder()
        b_pandora = ak.ArrayBuilder()
        b_cluster = ak.ArrayBuilder()   # ak array for MC cluster  1/2 if the cluster is from conversion/brems else 0, if both then 3
        b_event = ak.ArrayBuilder()     # ak array for event energy w/o neutrinos

        if idx < skip:
            print("Skipping event #", idx)
            continue
        if idx >= nread+skip:
            print("Reached maxread at idx ="+str(idx))
            break
        print("Reading",idx, "/", nread, "events...")

        # making dictionary for MC clusters
        dictClusterIdx = {}
        dictSimMc = {}
        # making dictionary for Pandora clusters
        dictPandora = {}    # calhit to pandora cluster id
        dictPandoraPrediction = {}    # pandora cluster id to pandora predictions and MC trith

        dictSkimmed = {}    # MCparticlesSkimmed
        skimmedIds = []     # skimmed Id
        dictDiJets ={}
        dictMergeSplit = {} # dictionary mcid to mcid (only modified in split (pi0->2gamma, brems))

        dictMCP = {}
        dictEvent = {}
        totalEnergy = 0.0
        di_quark_energy = 0.0

        dictTrackIds = {}

        nsimhit = 0
        colnames = event.getCollectionNames()

        for colname in SkimmedList:
            if not colname in colnames: continue
            col = event.getCollection(colname)
            for mcpSkimmed in col:
                skimmedIds.append(mcpSkimmed.id())
                if mcpSkimmed.getParents().size()==0:
                    dictSkimmed[mcpSkimmed.id()] = mcpSkimmed
                    continue
                mergedMCP = mcpSkimmed
                while mergedMCP.getParents().size() > 0:
                    # pi0 -> gamma gamma が近接していれば親(pi0)へマージ
                    if checkPi02gamma(mergedMCP):
                        parent = mergedMCP.getParents()[0]
                        photon1Pos = get_entry_point_from_simhits(parent.getDaughters()[0], event, SimHitNameList)
                        photon2Pos = get_entry_point_from_simhits(parent.getDaughters()[1], event, SimHitNameList)
                        # print(mergedMCP.id(), parent.getEnergy(), parent.getDaughters()[0].getEnergy(), parent.getDaughters()[1].getEnergy())
                        n_decayedpi0 = n_decayedpi0 + 1
                        if checkDistance(photon1Pos, photon2Pos):
                            mergedMCP = parent
                            n_decayedpi0_near = n_decayedpi0_near + 1
                            continue
                        break
                    # brems が近接していれば親(e-/e+)へマージ
                    if checkBrems(mergedMCP):
                        parent = mergedMCP.getParents()[0]
                        photonPos = get_entry_point_from_simhits(mergedMCP, event, SimHitNameList)
                        electronPos = get_calo_entry_from_mcp_via_track(parent, event, TrackList, RelTrackList)
                        if electronPos is None or not checkDistance(photonPos, electronPos):
                            n_brems = n_brems + 1
                            break
                        mergedMCP = parent
                        n_brems_near = n_brems_near + 1
                        continue
                    # どちらの条件にも当てはまらなければ終了
                    break

                dictSkimmed[mcpSkimmed.id()] = mergedMCP

        b_event.begin_list()

        for colname in MCParticleList:
            if not colname in colnames: continue
            col = event.getCollection(colname)
            for mcp in col:
                # print(mcp.getParents())
                if mcp.id() in skimmedIds:
                    dictSkimmed[mcp.id()] = mcp
                    # print(mcp.id(), mcp.id())
                else:
                    imcp = mcp
                    check_parentency = True
                    while imcp.getParents():
                        imcp = imcp.getParents()[0]
                        # print(mcp.id(), imcp.id())
                        if imcp.id() in skimmedIds:
                            dictSkimmed[mcp.id()] = imcp
                            check_parentency = False
                            break
                    if check_parentency:
                        print("  abnormal particle : ", mcp.id)

                if mcp.getDaughters().size()==0 and mcp.isCreatedInSimulation() and (mcp.getPDG() not in [12,14,16,18,-12,-14,-16,-18]): totalEnergy += mcp.getEnergy()
                if mcp.getParents().size()>0 and mcp.getParents()[0].getPDG()==25 and (mcp.getPDG() in [1,2,3,4,5,6, -1,-2,-3,-4,-5,-6]): 
                    print("  quark energy   ", mcp.getEnergy())
                    di_quark_energy += mcp.getEnergy()
                    dictDiJets[mcp.id()] = {"energy":mcp.getEnergy(), "pdg":mcp.getPDG(), "momentum":mcp.getMomentum()}
                    b_event.begin_list()
                    b_event.real(mcp.getEnergy())
                    b_event.real(mcp.getMomentum()[0])
                    b_event.real(mcp.getMomentum()[1])
                    b_event.real(mcp.getMomentum()[2])
                    b_event.end_list()
                
                dictMCP[mcp.id()] = returnDictMCP(mcp)
                # print(mcp.id(), dictMCP[mcp.id()])
                
        print("  jet energy     ", di_quark_energy)
        print("  non nu energy  ", totalEnergy)
        if len(dictDiJets)!=2: print("something is wrong at this event", idx)

        b_event.begin_list()
        b_event.real(totalEnergy)
        b_event.real(di_quark_energy)
        b_event.end_list()

        b_event.end_list()
        # print("totalEnergy", totalEnergy)

        for colname in SimHitNameList:
            if not colname in colnames: continue
            col = event.getCollection(colname)
            # print("col : ", col)
            for simhit in col:
                nsimhit += 1
                nmc = simhit.getNMCContributions()
                #if nmc > 1:
                #    print( "NMC: ", nmc, ", e1:", simhit.getEnergyCont(0), ", e2:", simhit.getEnergyCont(1))
                # currently only use the first term
                edepmax = 0
                for nmc in range(simhit.getNMCContributions()):
                    tmc = simhit.getParticleCont(nmc)
                    edep = simhit.getEnergyCont(nmc)
                    # print("  ", tmc.id())
                    if edepmax < edep:
                        edepmax = edep
                        # dictSimMc[simhit.id()] = tmc.id()
                        # mc = tmc
                        dictSimMc[simhit.id()] = dictSkimmed[tmc.id()].id()
                        mc = dictSkimmed[tmc.id()]
                    #print("tmcId : ",tmc.id(), "edep : ",edep)
                if edepmax == 0:
                    continue

                mcid = dictSimMc[simhit.id()]
                #### print("  ", simhit.id(), dictSimMc[simhit.id()])
                if not mcid in dictClusterIdx.keys():
                    #isDeltaRay(mc)
                    # bremsでenergyが下がっている場合は補正（同一MCの全labelで同じ値を使う）
                    # trackのbrems点での運動量を使うためevent等を渡す
                    corr_energy, corr_momentum = get_brems_corrected_energy_momentum(mc, event, TrackList, RelTrackList, dictSkimmed, colnames)
                    dictClusterIdx[mcid] = {"id":len(dictClusterIdx),"energy":corr_energy, "energy_mcp":mc.getEnergy(), "pdg":mc.getPDG(), "charge":mc.getCharge(), "mass":mc.getMass(), "momentum":corr_momentum,  "momentum_mcp":{mc.getMomentum()[0],mc.getMomentum()[1],mc.getMomentum()[2]}, "momentum_mcp_x":mc.getMomentum()[0], "momentum_mcp_y":mc.getMomentum()[1], "momentum_mcp_z":mc.getMomentum()[2], "momentum_mcp_at_end_x":mc.getMomentumAtEndpoint()[0],"momentum_mcp_at_end_y":mc.getMomentumAtEndpoint()[1],"momentum_mcp_at_end_z":mc.getMomentumAtEndpoint()[2], "status":mc.getSimulatorStatus(), "conversionBrems":dictMCP[mcid]}
                    # print("    ", mcid, "id", len(dictClusterIdx),"energy",mc.getEnergy(), "pdg",mc.getPDG(), "charge",mc.getCharge(), "mass",mc.getMass(), "momentum",mc.getMomentum()[0],mc.getMomentum()[1],mc.getMomentum()[2], "status",mc.getSimulatorStatus())
                    # print("                      ", "energy",mc.getEnergy(), "calc energy",math.sqrt(mc.getMass()**2+mc.getMomentum()[0]**2+mc.getMomentum()[1]**2+mc.getMomentum()[2]**2))

                    #if mc.getPDG() == 22 and mc.getParents().size() > 0 and mc.getParents()[0].getCharge() != 0:
                    #    pmc = mc.getParents()[0]
                    #    print("Photon with parent",pmc.getPDG(), "MC_ID", dictClusterIdx[mcid]["id"], "PMC_ID", pmc.id(), "photon momentum", mc.getMomentum()[0], mc.getMomentum()[1], mc.getMomentum()[2])

        for key, value in dictClusterIdx.items():
            if(value["energy"]!= value["energy_mcp"]): print(key, value)

        ncalhit = 0
        ncluster = 0
        for colname in PandoraList:
            if not colname in colnames: continue
            col = event.getCollection(colname)
            for hit in col:
                # print("Pandora ", hit, hit.id(), hit.getEnergy(), hit.getParticles(), hit.getClusters().size(), hit.getParticleIDs())
                pids = []
                for pid in hit.getParticleIDs():
                    pids.append([pid.getPDG(), pid.getLikelihood(), pid.getAlgorithmType()])
                for cluster in hit.getClusters():
                    ncluster += 1
                    clusterId = cluster.id()
                    # print("    cluster : ", cluster.id())

                    for relcolname in RelPandoraList:
                        if not relcolname in colnames: continue
                        relcol = event.getCollection(relcolname)
                        nav = LCRelationNavigator(relcol)
                        for rel,w in zip(nav.getRelatedToObjects(hit),nav.getRelatedToWeights(hit)):
                            if w > 0.5:
                                if w < 1: print (w)
                                # mcp = dictClusterIdx[rel.id()]
                                # print("    MCID",rel.id(),"assigned to",mcid)
                                # if mcp["energy"] != rel.getEnergy():
                                #     print("-------------------------------------------------------------", mcp["energy"], rel.getEnergy())
                                # # labels = dictClusterIdx[mcid]
                                # # hitid = hit.id()
                                break
                    if not clusterId in dictPandoraPrediction.keys():
                        dictPandoraPrediction[clusterId] = {"id":len(dictPandoraPrediction), "mcTruthClusterEnergy":rel.getEnergy(), "predictedPandoraPFOsEnergy":hit.getEnergy(), "predictedClusterEnergy":cluster.getEnergy(), "type":cluster.getType(), "pid":pids}
                        # print(hit.getEnergy(), cluster.getEnergy())

                    for calhit in cluster.getCalorimeterHits():
                        ncalhit += 1
                        # print("        calhit : ", calhit.id() ,clusterId, cluster.id())
                        calhitid = calhit.id()
                        if not calhitid in dictPandora.keys():
                            #isDeltaRay(mc)
                            # dictPandora[calhitid] = {"id":len(dictPandora), "clusterId":clusterId}#cluster.id()}#, "energy":hit.getEnergy(), "pdg":hit.getParticleIDs(), "charge":hit.getCharge(), "mass":hit.getMass(), "momentum":hit.getMomentum()}#, "status":hit.getSimulatorStatus()}
                            dictPandora[calhitid] = {"id":len(dictPandora), "clusterId":clusterId}#, "mcTruthClusterEnergy":rel.getEnergy(), "predictedClusterEnergy":cluster.getEnergy(), "type":cluster.getType(), "pid":pids}#, "energy":hit.getEnergy(), "pdg":hit.getParticleIDs(), "charge":hit.getCharge(), "mass":hit.getMass(), "momentum":hit.getMomentum()}#, "status":hit.getSimulatorStatus()}
                            # print(calhitid, "id",len(dictPandora), "clusterId", cluster.id(), clusterId)#, "energy":hit.getEnergy(), "pdg":hit.getParticleIDs(), "charge":hit.getCharge(), "mass":hit.getMass(), "momentum":hit.getMomentum()}#, "status":hit.getSimulatorStatus()}

        print("MC dictionary finished. #hits =", nsimhit, " # clusters =", len(dictClusterIdx))
        print("pandora dictionary finished. #hits =", ncalhit, " # clusters =", ncluster)#len(dictPandora))

        if nsimhit == 0:
            print("# simhit = 0; skipping event #", idx)
            continue
        if len(dictClusterIdx) == 0:
            print("# cluster = 0; skipping event #", idx)
            continue

        havehit=0

        # matching hits and collect features and labels
        for colname in CaloHitNameList:
            if not colname in colnames: continue
            col = event.getCollection(colname)
            for hit in col:
                #print(hit.id(),hit.getEnergy(), hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2], hit.getTime())

                if havehit == 0: # make entry only if having a hit
                    b_feat.begin_list() # list for event
                    b_label.begin_list()
                    b_pandora.begin_list()
                    havehit = 1

                with b_feat.list(): # for hit
                    b_feat.real(hit.getEnergy())
                    b_feat.real(hit.getPosition()[0])
                    b_feat.real(hit.getPosition()[1])
                    b_feat.real(hit.getPosition()[2])
                    b_feat.real(hit.getTime())
                    b_feat.real(0) # not a virtual hit representing a track
                    b_feat.real(0) # omega for track
                    b_feat.real(0) # px for track
                    b_feat.real(0) # py for track
                    b_feat.real(0) # pz for track
                    b_feat.real(0) # px_vertex for track
                    b_feat.real(0) # py_vertex for track
                    b_feat.real(0) # pz_vertex for track
                #print ("features finished.")

                b_label.begin_list()
                hitid=0
                for relcolname in RelHitList:
                    if not relcolname in colnames: continue
                    relcol = event.getCollection(relcolname)
                    nav = LCRelationNavigator(relcol)
                    for rel,w in zip(nav.getRelatedToObjects(hit),nav.getRelatedToWeights(hit)):
                        if w > 0.5:
                            if w < 1: print (w)
                            mcid = dictSimMc[rel.id()]
                            #print("MCID",rel.id(),"assigned to",mcid)
                            labels = dictClusterIdx[mcid]
                            hitid = hit.id()

                #if(labels["charge"] != 0):
                #    print("Hit for track: id", labels["id"], "pos", hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2],"time", hit.getTime())

                b_label.real(hitid)
                if(hitid != 0):
                    b_label.integer(labels["id"])
                    b_label.integer(labels["pdg"])
                    b_label.real(labels["charge"])
                    b_label.real(labels["mass"])
                    # print(labels["momentum"], labels["momentum"][0], labels["momentum"][1], labels["momentum"][2])
                    b_label.real(labels["momentum"][0])
                    b_label.real(labels["momentum"][1])
                    b_label.real(labels["momentum"][2])
                    b_label.real(labels["status"])
                else:
                    b_label.real(-1) #labels["id"])
                    b_label.real(0) #labels["pdg"])
                    b_label.real(0) #labels["charge"])
                    b_label.real(0) #labels["mass"])
                    b_label.real(0) #labels["momentum"][0])
                    b_label.real(0) #labels["momentum"][1])
                    b_label.real(0) #labels["momentum"][2])
                    b_label.real(0) #labels["status"])

                #print ("labels finished.")
                b_label.end_list()

                b_pandora.begin_list()
                b_pandora.real(hitid)
                if(hitid in dictPandora):
                    pandoras = dictPandora[hitid]
                    # print("  found hitid:", hitid, "clusterId", pandoras["clusterId"])
                else:
                    # print("  not found ", hitid)
                    pandoras = None
                if(pandoras):
                    pandoraPrediction = dictPandoraPrediction[pandoras["clusterId"]]
                    b_pandora.integer(pandoras["id"])
                    b_pandora.integer(pandoras["clusterId"])
                    b_pandora.real(pandoraPrediction["mcTruthClusterEnergy"])
                    # b_pandora.real(pandoraPrediction["predictedClusterEnergy"])
                    b_pandora.real(pandoraPrediction["predictedPandoraPFOsEnergy"])
                    b_pandora.real(pandoraPrediction["type"])
                    # b_pandora.real(pandoraPrediction["pid"])
                    b_pandora.real(0)
                    b_pandora.real(0)
                    b_pandora.real(0)
                    # b_pandora.real(pandoras["id"])
                    # b_pandora.real(pandoras["clusterId"])
                    # b_pandora.real(pandoras["energy"])
                    # b_pandora.real(0)
                    # b_pandora.real(pandoras["charge"])
                    # b_pandora.real(pandoras["mass"])
                    # b_pandora.real(0)
                    # b_pandora.real(0)
                    # if(hit.getPosition()[0]<0 and hit.getPosition()[0]>-1500 and hit.getPosition()[2]>1000 and hit.getPosition()[1]<-1000):
                    # print(hitid, pandoras["id"], pandoras["clusterId"]+1, hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2])
                    # if(pandoras["clusterId"]+1 == 6100):
                    #     print("      6100", hitid, pandoras["id"], pandoras["clusterId"]+1, hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2])
                    # print(pandoras["clusterId"]+1)
                else:
                    b_pandora.integer(-1) #labels["id"])
                    b_pandora.integer(-1) #labels["pdg"])
                    b_pandora.real(0) #labels["charge"])
                    b_pandora.real(0) #labels["mass"])
                    b_pandora.real(0) #labels["momentum"][0])
                    b_pandora.real(0) #labels["momentum"][1])
                    b_pandora.real(0) #labels["momentum"][2])
                    b_pandora.real(0) #labels["status"])
                    # print(0)
                    # print(hitid, -1, 0, hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2])
                b_pandora.end_list()

        if havehit == 0:
            print("No calohit in the event ", idx)
            continue

        ###################################################################
        # track -> MCParticle対応を一意化:
        # - 1つのMCParticleに複数trackが紐づく場合は、
        #   (merge時は decay 前MC を優先した上で) innermostHit半径が最大のtrackを採用
        # - 非採用trackは feat/label を作らない
        bestCandidateForTrack = {}
        groupedCandidatesForMC = {}

        for colname in TrackList:
            if not colname in colnames: continue
            col = event.getCollection(colname)
            for track in col:
                rin = track.getRadiusOfInnermostHit()
                if rin <= 0:
                    # fallback: track stateの参照点からinnermost(最小半径)を計算
                    rin = 1e30
                    for ts in track.getTrackStates():
                        ref = ts.getReferencePoint()
                        r = math.sqrt(ref[0]*ref[0] + ref[1]*ref[1] + ref[2]*ref[2])
                        if r < rin:
                            rin = r
                    if rin >= 1e30:
                        rin = -1.0
                if rin < 0:
                    continue

                bestRelForTrack = None
                for relcolname in RelTrackList:
                    if relcolname not in colnames: continue
                    relcol = event.getCollection(relcolname)
                    nav = LCRelationNavigator(relcol)
                    for rel, w in zip(nav.getRelatedToObjects(track), nav.getRelatedToWeights(track)):
                        if w <= 0.5: continue

                        # mergeしている場合は merge先(=decay前側)のMCIDを用いる
                        mapped_mcp = dictSkimmed[rel.id()] if rel.id() in dictSkimmed else rel
                        mapped_mcid = mapped_mcp.id()
                        if mapped_mcid not in dictClusterIdx:
                            continue

                        candidate = {
                            "track": track,
                            "mcid": mapped_mcid,
                            "weight": w,
                            "rin": rin,
                            "is_pre_decay": (rel.id() == mapped_mcid)
                        }
                        if bestRelForTrack is None:
                            bestRelForTrack = candidate
                        else:
                            # 同一trackに複数候補がある場合は、weight優先で代表を1つにする
                            if (candidate["weight"], candidate["is_pre_decay"], candidate["rin"]) > (
                                bestRelForTrack["weight"], bestRelForTrack["is_pre_decay"], bestRelForTrack["rin"]
                            ):
                                bestRelForTrack = candidate

                if bestRelForTrack is not None:
                    bestCandidateForTrack[track.id()] = bestRelForTrack
                    mcid = bestRelForTrack["mcid"]
                    if mcid not in groupedCandidatesForMC:
                        groupedCandidatesForMC[mcid] = []
                    groupedCandidatesForMC[mcid].append(bestRelForTrack)

        selectedTrackInfo = {}
        for mcid, candidates in groupedCandidatesForMC.items():
            print("MCID", mcid, "candidates", candidates, "track ids", [c["track"].id() for c in candidates])
            if len(candidates) == 1:
                winner = candidates[0]
            else:
                pre_decay_candidates = [c for c in candidates if c["is_pre_decay"]]
                pool = pre_decay_candidates if len(pre_decay_candidates) > 0 else candidates
                # innermostHitが最も外側(半径最大)のtrackを採用
                winner = max(pool, key=lambda c: (c["rin"], c["weight"], -c["track"].id()))
            selectedTrackInfo[winner["track"].id()] = {"mcid": mcid}
        print(selectedTrackInfo)

        for colname in TrackList:
            if not colname in colnames: continue
            col = event.getCollection(colname)
            for track in col:
                if track.id() not in selectedTrackInfo: continue
                nstates = track.getTrackStates().size()
                print("Track find. #states = ", nstates)
                if nstates > 0:
                    state = track.getTrackState(4) # at calorimeter
                    point = state.getReferencePoint()
                    #print("AtcaloPos: ", point[0], point[1], point[2])
                    p = -.001/state.getOmega() * 0.3 * 3.5 # Tesla
                    px = p * math.cos(state.getPhi())
                    py = p * math.sin(state.getPhi())
                    pz = p * state.getTanLambda()
                    #print("AtcaloMom: ", px, py, pz)

                    state_vertex = track.getTrackState(1) # at IP
                    point_vertex = state_vertex.getReferencePoint()
                    #print("AtcaloPos: ", point[0], point[1], point[2])
                    p_vertex = -.001/state_vertex.getOmega() * 0.3 * 3.5 # Tesla
                    px_vertex = p_vertex * math.cos(state_vertex.getPhi())
                    py_vertex = p_vertex * math.sin(state_vertex.getPhi())
                    pz_vertex = p_vertex * state_vertex.getTanLambda()
                    #print("AtcaloMom: ", px, py, pz)

                    with b_feat.list(): # record for features
                        b_feat.real(0.) # no edep in virtual track hit
                        b_feat.real(point[0])
                        b_feat.real(point[1])
                        b_feat.real(point[2])
                        b_feat.real(0.) # track have no time info
                        b_feat.real(1) # a virtual hit representing a track
                        # track-specific features
                        omega = track.getOmega()
                        b_feat.real((omega>0)-(omega<0)) # no sgn() in python (numpy has)
                        b_feat.real(px) # momentum at calo entry
                        b_feat.real(py)
                        b_feat.real(pz)
                        b_feat.real(px_vertex) # momentum at vertex or IP
                        b_feat.real(py_vertex)
                        b_feat.real(pz_vertex)

                    # b_label.begin_list()
                    # hitid=0
                    # 
                    # for relcolname in RelTrackList:
                    #     if not relcolname in colnames: continue
                    #     relcol = event.getCollection(relcolname)
                    #     nav = LCRelationNavigator(relcol)
                    #     for rel,w in zip(nav.getRelatedToObjects(track),nav.getRelatedToWeights(track)):
                    #         if w > 0.5:
                    #             if w < 1: print (w)
                    #             if rel.id() not in dictClusterIdx:
                    #                 print("MCID", rel.id(), "not found.")
                    #                 if rel.getEnergy()>1:
                    #                     print("     trackID : ", track.id(), " ,  energy : ", rel.getEnergy())
                    #                     if not rel.isStopped():
                    #                         notStopped = notStopped + 1
                    #                     if rel.isDecayedInTracker():
                    #                         decayedInTracker = decayedInTracker + 1
                    #                     abnormal_virtual_hit = abnormal_virtual_hit + 1
                    #                     if rel.getPDG() ==211 or rel.getPDG() == -211:
                    #                         npi = npi+1
                    #                     if rel.getPDG() == 321 or rel.getPDG() == -321:
                    #                         nk = nk+1
                    #                 continue
                    #             labels = dictClusterIdx[rel.id()]
                    #             hitid = -track.id()
                    # 
                    # b_label.real(hitid)

                    b_label.begin_list()
                    labels = None
                    hitid = 0
                    selected_mcid = selectedTrackInfo[track.id()]["mcid"]
                    if selected_mcid in dictClusterIdx:
                        labels = dictClusterIdx[selected_mcid]
                        hitid = -track.id()
                    else:
                        print("MCID", selected_mcid, "not found.")
                    b_label.real(hitid)

                    #print("Track pos", point[0], point[1], point[2], "momentum", px, py, pz)
                    if hitid != 0:
                    # if hitid == 0:
                        print("Track ", -hitid, " assigned to cluster ", labels["id"])

                    if(hitid != 0 and labels is not None):
                        b_label.real(labels["id"])
                        b_label.real(labels["pdg"])
                        b_label.real(labels["charge"])
                        b_label.real(labels["mass"])
                        b_label.real(labels["momentum"][0])
                        b_label.real(labels["momentum"][1])
                        b_label.real(labels["momentum"][2])
                        b_label.real(labels["status"])
                        if not labels["id"] in dictTrackIds.keys():
                            dictTrackIds[labels["id"]] = -hitid
                        else:
                            print("     clusterId ", labels["id"], " have multiple tracks ", dictTrackIds[labels["id"]], "  ", -hitid, "MC momentum ", labels["momentum"][0], labels["momentum"][1], labels["momentum"][2])
                    else:
                        b_label.real(-1) #labels["id"])
                        b_label.real(0) #labels["pdg"])
                        b_label.real(0) #labels["charge"])
                        b_label.real(0) #labels["mass"])
                        b_label.real(0) #labels["momentum"][0])
                        b_label.real(0) #labels["momentum"][1])
                        b_label.real(0) #labels["momentum"][2])
                        b_label.real(0) #labels["status"])
                        print("no matched MCParticle at track ", track.id())
                    b_label.end_list()

                    b_pandora.begin_list()
                    b_pandora.real(hitid)
                    b_pandora.real(0) #labels["id"])
                    b_pandora.real(-1) #labels["pdg"])
                    b_pandora.real(0) #labels["charge"])
                    b_pandora.real(0) #labels["mass"])
                    b_pandora.real(0) #labels["momentum"][0])
                    b_pandora.real(0) #labels["momentum"][1])
                    b_pandora.real(0) #labels["momentum"][2])
                    b_pandora.real(0) #labels["status"])
                    b_pandora.end_list()
                    # print(hitid, 0, 0, hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2])

        b_feat.end_list() # list for event
        b_label.end_list()
        b_pandora.end_list()

        b_cluster.begin_list()
        for key in dictClusterIdx.keys():
            _dictCluster = dictClusterIdx[key]
            b_cluster.begin_list()
            b_cluster.integer(_dictCluster["id"])
            b_cluster.integer(_dictCluster["conversionBrems"])
            b_cluster.end_list()
            # print(key, _dictCluster["id"], _dictCluster["conversionBrems"])
        b_cluster.end_list()
        # print(b_cluster)



        # snapshot
        feat = b_feat.snapshot()
        label = b_label.snapshot()
        pandora = b_pandora.snapshot()
        cluster = b_cluster.snapshot()
        event = b_event.snapshot()

        # np.set_printoptions(threshold=np.inf)
        # np_pandora = ak.to_numpy(pandora[:,0:100,2])
        # print(np_pandora)

        if array_exist:
            #print("concatenate")
            #print(ak_feat.type)
            ak_feat = ak.concatenate((ak_feat, feat), axis=0)
            ak_label = ak.concatenate((ak_label, label), axis=0)
            ak_pandora = ak.concatenate((ak_pandora, pandora), axis=0)
            ak_cluster = ak.concatenate((ak_cluster, cluster), axis=0)
            ak_event = ak.concatenate((ak_event, event), axis=0)
            #print(ak_feat.type)
        else:
            ak_feat = feat
            ak_label = label
            ak_pandora = pandora
            ak_cluster = cluster
            ak_event = event
            array_exist = 1


    # save to parquet
    #ak.to_parquet(ak_array,outfilename)

    #print(ak_feat.type)
    #print(ak_label.type)

    # to hdf5 file
    file = h5py.File(outfilename,"w")

    g_feat = file.create_group("feature")
    form, length, container = ak.to_buffers(ak_feat, container=g_feat)
    g_feat.attrs["form"] = form.to_json()
    g_feat.attrs["length"] = json.dumps(length)

    g_label = file.create_group("label")
    form, length, container = ak.to_buffers(ak_label, container=g_label)
    g_label.attrs["form"] = form.to_json()
    #g_label.attrs["form"] = ak.to_json(ak_label)
    g_label.attrs["length"] = json.dumps(length)

    g_pandora = file.create_group("pandora")
    form, length, container = ak.to_buffers(ak_pandora, container=g_pandora)
    g_pandora.attrs["form"] = form.to_json()
    g_pandora.attrs["length"] = json.dumps(length)

    g_cluster = file.create_group("cluster")
    form, length, container = ak.to_buffers(ak_cluster, container=g_cluster)
    g_cluster.attrs["form"] = form.to_json()
    g_cluster.attrs["length"] = json.dumps(length)

    g_event = file.create_group("event")
    form, length, container = ak.to_buffers(ak_event, container=g_event)
    g_event.attrs["form"] = form.to_json()
    g_event.attrs["length"] = json.dumps(length)

    print("notStopped : ", notStopped, "   decayedInTracker : ", decayedInTracker, "    abnormal_virtual_hit : ", abnormal_virtual_hit, "  pi ", npi, "   k ", nk)
    print("number of event : ", nread, "   total decayed pi0 : ", n_decayedpi0, "   combined decayed pi0 : ", n_decayedpi0_near)
    print("number of event : ", nread, "   total brems       : ", n_brems, "   combined brems       : ", n_brems_near)

    #ak0_array = ak.to_awkward0(ak_array)
    #awkward0.save(outfilename, ak0_array, mode="w")

    #s = pa.ipc.new_stream(outfilename,pa_batch.schema)
    #s.write_batch(pa_batch)

    #del s
    del reader
    return

#==========================================================================================
if __name__=='__main__':
    makeAk(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))


