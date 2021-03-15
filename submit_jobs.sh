dir=/mnt/home/kimjane7/MCRG/jobs

for job in "$dir"/*.sb
do
    sbatch "$job"
done

