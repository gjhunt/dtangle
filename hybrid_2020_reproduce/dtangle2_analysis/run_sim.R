L = list()
for(i in 1:nrow(grd)){
    params = grd[i,]
    L[[length(L)+1]] =
        tryCatch({
        run_sim(n_samples=50,mu=params$mu,
                                     sigma=params$sigma,
                  n_markers=params$n_markers,
                  seed=NULL,
                type=type,
                gen_sim=gen_sim)
        },error=function(e){
            print(e)
            "error"
        })
    print(Sys.time())
    print(i)
}

print("saving RDS...")
saveRDS(list(L=L,grd=grd,type=type),fname)
