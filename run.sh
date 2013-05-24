NCPU=4
WAVE3D_K=64
GEOM=SubmarineJ
BASE=$GEOM.wrl_${WAVE3D_K}_20_$NCPU

mpirun -np $NCPU ./tt \
-posfile $BASE/pos \
-denfile $BASE/den \
-geomprtn $BASE/geomprtn \
-valfile $BASE/val \
-chkfile $BASE/chk \
-knl 0 \
-mlib3d_NPQ 4 \
-mlib3d_ldname \
helm3d_ld_1.bin \
-mlib3d_hdname \
helm3d_hd_1_4.bin \
-wave3d_ACCU 1 \
-wave3d_NPQ 4 \
-wave3d_K $WAVE3D_K \
-wave3d_ptsmax 80 \
-wave3d_maxlevel 12