using Distributed, ClusterManagers
data_root = "/local/grahams"

function init_fast_partition(ntasks=24; partition="fast",
                                        time="1-00:00:00",
                                        mem_per_cpu=5000,
                                        mail_user="remote.compute.results@gmail.com",
                                        mail_type="END",
                                        kwargs...)
    addprocs_slurm(ntasks; partition=partition,
                           time=time,
                           mem_per_cpu=mem_per_cpu,
                           mail_user=mail_user,
                           mail_type=mail_type,
                           kwargs...)
end
    
function init_partition(ntasks=1; partition="debug",
                                        time="08:00:00",
                                        mem_per_cpu=5000,
                                        mail_user="remote.compute.results@gmail.com",
                                        mail_type="END",
                                        kwargs...)
    addprocs_slurm(ntasks; partition=partition,
                           time=time,
                           mem_per_cpu=mem_per_cpu,
                           mail_user=mail_user,
                           mail_type=mail_type,
                           kwargs...)
end
    
