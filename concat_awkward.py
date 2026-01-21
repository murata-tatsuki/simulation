import json
import h5py
import awkward as ak
import numpy as np
import sys
import load_awkward as la

def load_awkward(filename):
    file = h5py.File(filename,"r")
    group = file["awkward"]
    
    form = ak.forms.from_json(group.attrs["form"])
    length = json.loads(group.attrs["length"])
    container = {k: np.asarray(v) for k, v in group.items()}

    return ak.from_buffers(form, length, container)    

def load_awkward2(filename):
    file = h5py.File(filename,"r")
    feat = file["feature"]
    
    form = ak.forms.from_json(feat.attrs["form"])
    length = json.loads(feat.attrs["length"])
    container = {k: np.asarray(v) for k, v in feat.items()}

    ak_feat = ak.from_buffers(form, length, container)    

    label = file["label"]
    
    form = ak.forms.from_json(label.attrs["form"])
    length = json.loads(label.attrs["length"])
    container = {k: np.asarray(v) for k, v in label.items()}

    ak_label = ak.from_buffers(form, length, container)    

    return ak_feat, ak_label

def load_awkward2_pandora(filename):
    file = h5py.File(filename,"r")
    feat = file["feature"]
    
    form = ak.forms.from_json(feat.attrs["form"])
    length = json.loads(feat.attrs["length"])
    container = {k: np.asarray(v) for k, v in feat.items()}
    ak_feat = ak.from_buffers(form, length, container)    


    label = file["label"]
    
    form = ak.forms.from_json(label.attrs["form"])
    length = json.loads(label.attrs["length"])
    container = {k: np.asarray(v) for k, v in label.items()}
    ak_label = ak.from_buffers(form, length, container)    


    label = file["pandora"]
    
    form = ak.forms.from_json(label.attrs["form"])
    length = json.loads(label.attrs["length"])
    container = {k: np.asarray(v) for k, v in label.items()}
    ak_pandora = ak.from_buffers(form, length, container)    

    return ak_feat, ak_label, ak_pandora

def load_awkwards(filenames):
    for i, file in enumerate(filenames):
        feat, label = load_awkward2(file)
        if i==0:
            ak_feats = feat
            ak_labels = label
        else:
            ak_feats = ak.concatenate((ak_feats, feat), axis=0)
            ak_labels = ak.concatenate((ak_labels, label), axis=0)

        #print (ak.num(ak_feats,axis=0), ak.num(ak_labels,axis=0))

    return ak_feats, ak_labels

def load_awkwards_pandora(filenames, pandora=False):
    for i, file in enumerate(filenames):
        if(pandora):
            feat, label, pand = load_awkward2_pandora(file)
            if i==0:
                ak_feats = feat
                ak_labels = label
                ak_pandora = pand

            else:
                ak_feats = ak.concatenate((ak_feats, feat), axis=0)
                ak_labels = ak.concatenate((ak_labels, label), axis=0)
                ak_pandora = ak.concatenate((ak_pandora, pand), axis=0)

        #print (ak.num(ak_feats,axis=0), ak.num(ak_labels,axis=0))


        else: 
            feat, label = load_awkward2(file)
            if i==0:
                ak_feats = feat
                ak_labels = label
            else:
                ak_feats = ak.concatenate((ak_feats, feat), axis=0)
                ak_labels = ak.concatenate((ak_labels, label), axis=0)

        #print (ak.num(ak_feats,axis=0), ak.num(ak_labels,axis=0))

    if(pandora):
        return ak_feats, ak_labels, ak_pandora
    else:
        return ak_feats, ak_labels

def save_awkward(ak_feat, ak_label, outfilename):
    # to hdf5 file
    file = h5py.File(outfilename,"w")
    g_feat = file.create_group("feature")

    form, length, container = ak.to_buffers(ak_feat, container=g_feat)
    g_feat.attrs["form"] = form.to_json()
    g_feat.attrs["length"] = json.dumps(length)

    g_label = file.create_group("label")

    form, length, container = ak.to_buffers(ak_label, container=g_label)
    g_label.attrs["form"] = form.to_json()
    g_label.attrs["length"] = json.dumps(length)

def save_awkward_pandora(ak_feat, ak_label, ak_pandora, outfilename, pandora=False):
    # to hdf5 file
    file = h5py.File(outfilename,"w")
    g_feat = file.create_group("feature")

    form, length, container = ak.to_buffers(ak_feat, container=g_feat)
    g_feat.attrs["form"] = form.to_json()
    g_feat.attrs["length"] = json.dumps(length)

    g_label = file.create_group("label")

    form, length, container = ak.to_buffers(ak_label, container=g_label)
    g_label.attrs["form"] = form.to_json()
    g_label.attrs["length"] = json.dumps(length)

    g_pandora = file.create_group("pandora")

    form, length, container = ak.to_buffers(ak_pandora, container=g_pandora)
    g_pandora.attrs["form"] = form.to_json()
    g_pandora.attrs["length"] = json.dumps(length)



# main
filelist = sys.argv[1]
outfilename = sys.argv[2]

f = open(filelist,"r")
files = []

for x in f:
    files.append(x.rstrip("\n"))
f.close()

# bool_pandora=True
# if(bool_pandora):
#     ak_feat, ak_label, ak_pandora = load_awkwards_pandora(files, bool_pandora)
#     save_awkward_pandora(ak_feat, ak_label, ak_pandora, outfilename, bool_pandora)
# else:
#     ak_feat, ak_label = load_awkwards(files)
#     save_awkward(ak_feat, ak_label, outfilename)

ak_feats, ak_labels, ak_preds, ak_energys, ak_pandoras, ak_clusters, ak_events = la.load_awkwards_all(files)
la.save_awkward(outfilename, ak_feats, ak_labels, ak_pred=ak_preds, ak_energy=ak_energys, ak_x=None, ak_y=None, ak_pandora=ak_pandoras, ak_eventEnergy = ak_events, ak_bremsConversion = ak_clusters)

