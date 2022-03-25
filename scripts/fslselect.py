#!/usr/bin/env python

import nibabel as nib
import numpy as np
import sys
import re
import argparse

def argument_parser_fslselect(argv):
    parser=argparse.ArgumentParser(description='Select a set of ROI labels from a volume')
    parser.add_argument('-seq',action='store_true',help='Output ROIs are labeled sequentially in the order provided')
    parser.add_argument('-bin',action='store_true',help='Output is binarized')
    parser.add_argument('invol',action='store')
    parser.add_argument('outvol',action='store')
    parser.add_argument('roivals',action='store',nargs='*',help="val1 val2 val3=newval val4-val7 val8-val10=20 ...")
    
    return parser.parse_args(argv)

def parselist(liststring):
    outlist=[]
    for p in liststring.split(","):
        if not p:
            continue
        if "-" in p:
            pp=p.split("-")
            outlist+=np.arange(int(pp[0]),int(pp[-1])+1).tolist()
        else:
            outlist+=[int(p)]
    return outlist

def main(argv):
    args=argument_parser_fslselect(argv)
    
    isseq=args.seq
    isbin=args.bin
    involfile=args.invol
    outvolfile=args.outvol
    roistrings=args.roivals
    
    V=nib.load(involfile)
    Vdata=V.get_fdata()
    
    seqlist=[]
    seqlist_newval=[]
    maxval=0
    
    #True if seq OR if we ever assigned a new value using "="
    isnewval=isseq
    
    for v in roistrings:
      vve=v.split("=")
      newval=None
      if len(vve)==2:
          newval=int(vve[1])
          isnewval=True
      #vl=list(filter(None,re.split("[,+]",vve[0])))
      
      
      vv=parselist(vve[0])
      seqlist+=vv
      
      if newval is None:
          if isseq:
              seqlist+=np.arange(maxval+1,maxval+len(vv)).tolist()
          else:
              seqlist_newval+=vv
      else:
          seqlist_newval+=[newval]*len(vv)
          
      maxval=max(seqlist_newval)
      
    Vmask=np.isin(Vdata,seqlist)
    
    if isseq or isnewval:
        uvals,uinv=np.unique(Vdata[Vmask],return_inverse=True)
        #u2seq=np.isin(seqlist,uvals)
        #u2seqidx=np.where(u2seq)[0]
        
        u2seqidx=np.array([np.where(seqlist==x)[0][0] for x in uvals])
        Vnewdata=np.zeros(Vdata.shape)
        #Vnewdata[Vmask]=u2seqidx[uinv]+1
        
        Vnewdata[Vmask]=np.array(seqlist_newval)[u2seqidx[uinv]]
    else:
        Vnewdata=Vdata*Vmask
    
    if isbin:
        Vnewdata=Vnewdata>0
    
    if isinstance(V,nib.nifti1.Nifti1Image):
        Vnew=nib.Nifti1Image(Vnewdata.astype(V.get_data_dtype()),affine=V.affine, header=V.header)
    elif isinstance(V,nib.freesurfer.mghformat.MGHImage):
        Vnew=nib.MGHImage(Vnewdata.astype(V.get_data_dtype()),affine=V.affine, header=V.header)
    nib.save(Vnew,outvolfile)

if __name__ == "__main__":
    main(sys.argv[1:])
    
