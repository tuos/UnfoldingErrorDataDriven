
cd /home/tuos/workingpA/CMSSW_5_3_20/src
#cmsenv
eval `scramv1 runtime -sh`

cd /home/tuos/workingpA/v1/Unfold


export coll="PPb"
export MIN=-1.0 
export MAX=1.0 

  echo "root -l -b -q pAUnfoldSmearedMatrix.C++"

root -b > run.log <<EOF
.L pAUnfoldSmearedMatrix.C+
pAUnfoldSmearedMatrix(3, 0, 0, 0, 1, "", 20, 0, 4, 0, 0)
.q
EOF

export MIN=-2.5
export MAX=-2.0

  echo "root -l -b -q pAUnfoldSmearedMatrix.C++"

root -b > run.log <<EOF
.L pAUnfoldSmearedMatrix.C+
pAUnfoldSmearedMatrix(3, 0, 0, 0, 1, "", 20, 0, 4, 0, 0)
.q
EOF

export MIN=-2.0
export MAX=-1.5

  echo "root -l -b -q pAUnfoldSmearedMatrix.C++"

root -b > run.log <<EOF
.L pAUnfoldSmearedMatrix.C+
pAUnfoldSmearedMatrix(3, 0, 0, 0, 1, "", 20, 0, 4, 0, 0)
.q
EOF

export MIN=-1.5
export MAX=-1.0

  echo "root -l -b -q pAUnfoldSmearedMatrix.C++"

root -b > run.log <<EOF
.L pAUnfoldSmearedMatrix.C+
pAUnfoldSmearedMatrix(3, 0, 0, 0, 1, "", 20, 0, 4, 0, 0)
.q
EOF

export MIN=-1.0
export MAX=-0.5

  echo "root -l -b -q pAUnfoldSmearedMatrix.C++"

root -b > run.log <<EOF
.L pAUnfoldSmearedMatrix.C+
pAUnfoldSmearedMatrix(3, 0, 0, 0, 1, "", 20, 0, 4, 0, 0)
.q
EOF

export MIN=-0.5
export MAX=0.5

  echo "root -l -b -q pAUnfoldSmearedMatrix.C++"

root -b > run.log <<EOF
.L pAUnfoldSmearedMatrix.C+
pAUnfoldSmearedMatrix(3, 0, 0, 0, 1, "", 20, 0, 4, 0, 0)
.q
EOF

export MIN=0.5
export MAX=1.0

  echo "root -l -b -q pAUnfoldSmearedMatrix.C++"

root -b > run.log <<EOF
.L pAUnfoldSmearedMatrix.C+
pAUnfoldSmearedMatrix(3, 0, 0, 0, 1, "", 20, 0, 4, 0, 0)
.q
EOF

export MIN=1.0
export MAX=1.5

  echo "root -l -b -q pAUnfoldSmearedMatrix.C++"

root -b > run.log <<EOF
.L pAUnfoldSmearedMatrix.C+
pAUnfoldSmearedMatrix(3, 0, 0, 0, 1, "", 20, 0, 4, 0, 0)
.q
EOF



#echo "Copying output files to " $destination
