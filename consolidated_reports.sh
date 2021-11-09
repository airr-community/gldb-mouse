for strain in CAST_EiJ LEWES_EiJ MSM_MsJ NOD_ShiLtJ PWD_PhJ
do
cd $strain
python "../python/consolidate_mouse.py" "db/mouse_"$strain"_db.csv" ../imgt/musculus_imgt_ungapped.fasta "$strain"_summary.csv -source "watson_et_al/imcb_"$strain".fasta"
cd ..
done