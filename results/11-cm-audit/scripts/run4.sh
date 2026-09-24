cd /mnt/raid10/sim-work/tps62933/cm_audit
TB=0.00143
for V in fullcu topcut dualcut; do
  for TG in air tg; do
    if [ $TG = tg ]; then OPT="--topgnd"; else OPT=""; fi
    TAG=_${V}_${TG}
    python3 gen_geo4.py --var=$V --tb=$TB $OPT --tag=$TG >/dev/null
    gmsh -3 geo4${TAG}.geo -o m4${TAG}.msh -v 2 >/dev/null 2>&1 || { echo "GMSH FAIL $TAG"; continue; }
    ElmerGrid 14 2 m4${TAG}.msh -autoclean -out mesh4${TAG} >/dev/null 2>&1 || { echo "GRID FAIL $TAG"; continue; }
    echo "MESH $TAG names: $(tr "\n" " " < mesh4${TAG}/mesh.names)"
    for ER in 1.0 4.3 4.6; do
      RD=run4${TAG}_e${ER}
      rm -rf $RD; mkdir -p $RD
      python3 make_sif4.py mesh4${TAG} $ER $TB >/dev/null
      cp case4_mesh4${TAG}.sif $RD/case.sif
      cp -r mesh4${TAG} $RD/
      ( cd $RD && timeout 900 ElmerSolver case.sif > run.log 2>&1 )
      C=$(awk "NR==1{print \$4}" $RD/scalars.mesh4${TAG}.dat 2>/dev/null)
      echo "RESULT $TAG epsr=$ER C=$C F"
    done
  done
done
echo GEO4DONE
