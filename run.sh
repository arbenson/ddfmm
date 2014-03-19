# Distributed Directional Fast Multipole Method
#   Copyright (C) 2013 Austin Benson, Lexing Ying, and Jack Poulson
#
# This file is part of DDFMM.
#
#    DDFMM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    DDFMM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with DDFMM.  If not, see <http://www.gnu.org/licenses/>.

# Parameters
NCPU=2
WAVE3D_K=32
NPW=10
GEOM=sphere
BASE=${GEOM}.wrl_${WAVE3D_K}_${NPW}_${NCPU}
ACCU=2

# Run program
mpirun \
-np ${NCPU} ./tt \
-posfile ${BASE}/pos \
-denfile ${BASE}/den \
-geomprtn ${BASE}/geomprtn \
-geomfile ${GEOM}.wrl \
-valfile ${BASE}/val \
-chkfile ${BASE}/chk \
-knl 0 \
-mlib3d_NPQ 4 \
-mlib3d_ldname helm3d_ld_${ACCU}.bin \
-mlib3d_hdname helm3d_hd_${ACCU}_4.bin \
-wave3d_ACCU ${ACCU} \
-wave3d_NPQ 4 \
-wave3d_K ${WAVE3D_K} \
-wave3d_NPW ${NPW} \
-wave3d_ptsmax 80 \
-wave3d_maxlevel 12
