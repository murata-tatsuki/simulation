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

    for idx,event in enumerate(reader):
        b_feat = ak.ArrayBuilder()
        b_label = ak.ArrayBuilder()
        b_pandora = ak.ArrayBuilder()

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

        nsimhit = 0
        colnames = event.getCollectionNames()

        for colname in SkimmedList:
            if not colname in colnames: continue
            col = event.getCollection(colname)
            for mcpSkimmed in col:
                skimmedIds.append(mcpSkimmed.id())

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
                    dictClusterIdx[mcid] = {"id":len(dictClusterIdx),"energy":mc.getEnergy(), "pdg":mc.getPDG(), "charge":mc.getCharge(), "mass":mc.getMass(), "momentum":mc.getMomentum(), "status":mc.getSimulatorStatus()}
                    #### print("    ", mcid, "id", len(dictClusterIdx),"energy",mc.getEnergy(), "pdg",mc.getPDG(), "charge",mc.getCharge(), "mass",mc.getMass(), "momentum",mc.getMomentum(), "status",mc.getSimulatorStatus())

                    #if mc.getPDG() == 22 and mc.getParents().size() > 0 and mc.getParents()[0].getCharge() != 0:
                    #    pmc = mc.getParents()[0]
                    #    print("Photon with parent",pmc.getPDG(), "MC_ID", dictClusterIdx[mcid]["id"], "PMC_ID", pmc.id(), "photon momentum", mc.getMomentum()[0], mc.getMomentum()[1], mc.getMomentum()[2])

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


        for colname in TrackList:
            if not colname in colnames: continue
            col = event.getCollection(colname)
            for track in col:
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

                    b_label.begin_list()
                    hitid=0

                    for relcolname in RelTrackList:
                        if not relcolname in colnames: continue
                        relcol = event.getCollection(relcolname)
                        nav = LCRelationNavigator(relcol)
                        for rel,w in zip(nav.getRelatedToObjects(track),nav.getRelatedToWeights(track)):
                            if w > 0.5:
                                if w < 1: print (w)
                                if rel.id() not in dictClusterIdx:
                                    print("MCID", rel.id(), "not found.")
                                    if rel.getEnergy()>1:
                                        print("     trackID : ", track.id(), " ,  energy : ", rel.getEnergy())
                                        if not rel.isStopped():
                                            notStopped = notStopped + 1
                                        if rel.isDecayedInTracker():
                                            decayedInTracker = decayedInTracker + 1 
                                        abnormal_virtual_hit = abnormal_virtual_hit + 1
                                        if rel.getPDG() ==211 or rel.getPDG() == -211:
                                            npi = npi+1
                                        if rel.getPDG() == 321 or rel.getPDG() == -321:
                                            nk = nk+1
                                    continue
                                labels = dictClusterIdx[rel.id()]
                                hitid = -track.id()

                    b_label.real(hitid)

                    #print("Track pos", point[0], point[1], point[2], "momentum", px, py, pz) 
                    # if hitid != 0:
                    if hitid == 0:
                        print("Track ", -hitid, " assigned to cluster ", labels["id"])

                    if(hitid != 0):
                        b_label.real(labels["id"])
                        b_label.real(labels["pdg"])
                        b_label.real(labels["charge"])
                        b_label.real(labels["mass"])
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
        

            
        # snapshot
        feat = b_feat.snapshot()
        label = b_label.snapshot()
        pandora = b_pandora.snapshot()

        # np.set_printoptions(threshold=np.inf)
        # np_pandora = ak.to_numpy(pandora[:,0:100,2])
        # print(np_pandora)

        if array_exist:
            #print("concatenate")
            #print(ak_feat.type)
            ak_feat = ak.concatenate((ak_feat, feat), axis=0)
            ak_label = ak.concatenate((ak_label, label), axis=0)
            ak_pandora = ak.concatenate((ak_pandora, pandora), axis=0)
            #print(ak_feat.type)
        else:
            ak_feat = feat
            ak_label = label
            ak_pandora = pandora
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

    print("notStopped : ", notStopped, "   decayedInTracker : ", decayedInTracker, "    abnormal_virtual_hit : ", abnormal_virtual_hit, "  pi ", npi, "   k ", nk)


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
    
    
