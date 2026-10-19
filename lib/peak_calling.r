library(patchwork)

calc_window_count = function(df, window_center, window_size = 30){
    df %>%
        filter(TTS <= window_center + window_size / 2,
               TTS >= window_center - window_size / 2) %>% 
        nrow()
                                        # group_by(TTS) %>% 
                                        # summarise(n = n(), .groups = "drop")
}

calc_count_in_window = function(df, window_size = 30){
    
    list_window_count = list()
    for(i in unique(df$TTS)){
        list_window_count[[length(list_window_count)+1]] = 
            data.frame(TTS = i,
                       window_count = calc_window_count(df,
                                                        window_center = i, 
                                                        window_size = window_size))
    }
    
    df %>% 
        dplyr::left_join(
                   bind_rows(list_window_count),
                   by = "TTS"
               )
}


remove_noise = function(df){
    
    df_grouped = df %>% 
        group_by(TTS) %>% 
        summarise(n = n(), 
                  .groups = "drop")
    
    tmp_noise_level = 
        max(2, mean(df_grouped$n))
    
    df_grouped_filtered = 
        df_grouped %>% 
        mutate(n = floor(n - tmp_noise_level)) %>% 
        filter(n < 0)
                                        #return(df_grouped_filtered)
    tmp_list = list()
    for(i in unique(df_grouped_filtered$TTS)){
                                        # print(i)
        tmp_num = df_grouped_filtered %>% 
            filter(TTS == i) %>% 
            .$n
        
                                        # print(tmp_num)
        
        tmp_list[[length(tmp_list)+1]] = df %>% 
            filter(TTS == i) %>% 
            head(tmp_num)
    }
    
    return(bind_rows(tmp_list))
}

robust_calc_diff = function(x){
    tmp_diff = diff(x)
    mean_diff = mean(abs(tmp_diff))
    return(sapply(tmp_diff,
                  function(x){
                      if(abs(x) < mean_diff * 1.001){
                          return(0)
                      }else{
                          return(x)
                      }
                  }))
}

find_maxima = function(x, kde_bw = 30, kde_min_gradient_length = 0){
    
    stopifnot(length(x)>0)
    
    if(length(unique(x)) == 1){
        return(list(
            local_maxima = unique(x),
            global_maxima = unique(x)
        ))
    }
    
    if(max(x) - min(x) < 30){
        occurrence = sapply(unique(x), function(y){sum(x == y)})
        return(list(
            local_maxima = unique(x)[occurrence == max(occurrence)][1],
            global_maxima = unique(x)[occurrence == max(occurrence)][1]
        ))
    }
    
    tmp_density = density(x, kernel = "gaussian", bw = kde_bw, 
                          adjust = 1, 
                          n = max(x) - min(x)
                                        # n = 100
                          )
    
    local_optima_pos = cumsum(rle(sign(robust_calc_diff(tmp_density$y)))$lengths)
                                        # print("tmp_density$y")
                                        # print(tmp_density$y)
                                        # print("robust calc diff")
                                        # print(robust_calc_diff(tmp_density$y))
    ## find local maximum
    tmp_optima_i = c()
    tmp_values = rle(sign(robust_calc_diff(tmp_density$y)))$values
    tmp_lengths = rle(sign(robust_calc_diff(tmp_density$y)))$lengths
    
                                        # print(paste0("tmp_density$y: ", paste0(tmp_density$y, collapse = ", ")))
                                        # print(paste0("density rle values: ", paste0(tmp_values, collapse = ", ")))
                                        # print(paste0("density rle length: ", paste0(tmp_lengths, collapse = ", ")))
                                        # print(tmp_lengths)
                                        # print(length(tmp_values))
                                        # print(paste0("bw.JS ", tmp_bw))
    
    ## todo kde_min_gradient_length shall equal to the average length of 
    for(i in seq(length(tmp_values))){
                                        # if(i == 1 | i == length(tmp_values)){
                                        #   tmp_optima_i = c(tmp_optima_i, i)
                                        #   next
                                        # }
        if(tmp_values[i] == 1 &&
           tmp_lengths[i] > kde_min_gradient_length){
            if(i+1 > length(tmp_values)){
                next
            }
            for(j in seq(i+1, length(tmp_values))){
                                        # print(paste0(i, ", ",j ))
                if(tmp_values[j] == 0){
                    
                }else if(tmp_values[j] == -1 &&
                         tmp_lengths[j] > kde_min_gradient_length){
                    tmp_optima_i = c(tmp_optima_i, i)
                    break
                }else{
                    next
                }
            }
        }
    }
    
                                        # print(local_optima_pos)
                                        # print(tmp_optima_i)
    
    local_maxima = local_optima_pos[tmp_optima_i]
    local_maxima_x = tmp_density$x[local_maxima]
    local_maxima_y = tmp_density$y[local_maxima]
    local_maxima_x_filtered = local_maxima_x[tmp_density$y[local_maxima] == max(tmp_density$y[local_maxima])]
    
    return(list(local_maxima = round(local_maxima_x),
                global_maxima = round(local_maxima_x_filtered)))
}

find_peaks = function(df_reads, pA = "+", strand = "+", 
                      kde_bw = 30, 
                      kde_min_gradient_length = 10,
                      window_size = 5, 
                      window_count_filter = 1, 
                      pseudocount = 1){
    
                                        # stopifnot(nrow(df_reads) > 0)
    if(nrow(df_reads) == 0){
        return(data.frame())
    }
                                        # print(df_reads)
    df_reads_labeled = df_reads %>% 
        mutate(TTS = ifelse(V6 == "+", V8, V7),
               pA_boolean = ifelse(V5 > 15, "+", "-")) %>% 
        filter(V6 == strand,
               pA_boolean == pA) 
    
                                        # stopifnot(nrow(df_reads_labeled) > 0)
    if(nrow(df_reads_labeled) == 0){
        return(data.frame())
    }
    
    if(window_count_filter < 0){
        df_reads_filered =  df_reads_labeled
    }else{
        df_reads_filered = df_reads_labeled %>% 
            calc_count_in_window(df = .,
                                 window_size = window_size) %>% 
            filter(window_count > window_count_filter)
    }
    
    real.max.df <- df_reads_filered %>%
        group_by(TTS) %>%
        summarize(counts = n())
                                        # return(df_reads_filered)
                                        # filter(window_count > window_count_filter)
    
                                        # tmp_TTS_add_pseudo = c(df_reads_filered$TTS, rep())
    if(nrow(df_reads_filered) == 0){
        return(data.frame())
    }
                                        # stopifnot(nrow(df_reads_filered) > 0)
                                        # print(df_reads_filered)
    
    res_maxima = find_maxima(df_reads_filered$TTS, kde_bw = kde_bw, kde_min_gradient_length = kde_min_gradient_length)
    
    if(length(res_maxima$local_maxima)==0){
        warning("Cannot find a peak using the current parameter.\nPlease consider loose the window count filter, kde_min_gradient_length, or decrease kde_bw.")
        ## print(df_reads_filered)
        return(
            data.frame(
                pA = pA,
                strand = strand,
                local_maxima = NA,
                local_maxima.adj = NA, 
                local_maxima_window_counts = NA,
                global_maxima = NA,
                global_maxima_window_count = NA,
                total_count = nrow(df_reads_labeled)
            )
        )
    }
    local_maxima = res_maxima$local_maxima
    global_maxima = res_maxima$global_maxima

    global_count = calc_window_count(df_reads_labeled,
                                     window_center = global_maxima, 
                                     window_size = window_size)

    local_counts = sapply(local_maxima,
                          function(x){
                              calc_window_count(df_reads_labeled,
                                                window_center = x, 
                                                window_size = window_size)
                          })
    local_maxima.adj <- lapply(local_maxima,
                               function(x){ a <- real.max.df %>%
                                                dplyr::filter(TTS < x + (window_size/2) & TTS > x - (window_size/2))  %>%
                                                .[which.max(.$counts),] %>%
                                                .$TTS
                                                if(length(a) == 0) {
                                                    return(x)
                                                }else{
                                                    return(a)
                                                }
                               }
                               )

    return(data.frame(
        pA = pA,
        strand = strand,
        local_maxima = paste0(local_maxima, collapse = ","),
        local_maxima.adj = paste0(local_maxima.adj, collapse = ","),
        local_maxima_window_counts = paste0(local_counts, collapse = ","),
        global_maxima = global_maxima,
        global_maxima_window_count = global_count,
        total_count = nrow(df_reads_labeled)
    ))
                                        # print(paste0("local maxima: ", paste0(local_maxima, collapse = ", ")))
                                        # print(paste0("global maxima: ", global_maxima))
                                        # 
                                        # g1 = df_reads_labeled %>% 
                                        #   ggplot(aes(x = TTS))+
                                        #   geom_histogram(bins = 100)+
                                        #   geom_vline(xintercept = local_maxima[local_maxima != global_maxima], color = "cyan")+
                                        #   geom_vline(xintercept = global_maxima, color = "hotpink")+
                                        #   labs(title = paste0("strand ",strand))+
                                        #   # scale_y_log10()+
                                        #   # facet_wrap(~ V14)+
                                        #   theme_bw()
                                        # 
                                        # 
                                        # g2 = df_reads_labeled %>% 
                                        #   ggplot(aes(x = TTS))+
                                        #   geom_histogram(bins = 100)+
                                        #   geom_vline(xintercept = local_maxima[local_maxima != global_maxima], color = "cyan")+
                                        #   geom_vline(xintercept = global_maxima, color = "hotpink")+
                                        #   # scale_y_log10()+
                                        #   # facet_wrap(~ V14)+
                                        #   theme_bw()
                                        # 
                                        # if(length(global_maxima) > 0){
                                        #   g2 = g2+
                                        #     lims(x = c(global_maxima - 1000, global_maxima+1000))
                                        # }
                                        # return(g1 / g2)
}


show_peaks = function(df_reads, pA = "+", strand = "+", 
                      kde_bw = 30, 
                      kde_min_gradient_length = 10,
                      window_size = 5, 
                      window_count_filter = 1, 
                      pseudocount = 1){
    
                                        # stopifnot(nrow(df_reads) > 0)
    if(nrow(df_reads) == 0){
        return(NULL)
    }
                                        # print(df_reads)
    df_reads_labeled = df_reads %>% 
        mutate(TTS = ifelse(V6 == "+", V8, V7),
               pA_boolean = ifelse(V5 > 15, "+", "-")) %>% 
        filter(V6 == strand,
               pA_boolean == pA) 
    
                                        # stopifnot(nrow(df_reads_labeled) > 0)
    if(nrow(df_reads_labeled) == 0){
        return(NULL)
    }
    
    if(window_count_filter < 0){
        df_reads_filered =  df_reads_labeled
    }else{
        df_reads_filered = df_reads_labeled %>% 
            calc_count_in_window(df = .,
                                 window_size = window_size) %>% 
            filter(window_count > window_count_filter)
    }
    
    
                                        # return(df_reads_filered)
                                        # filter(window_count > window_count_filter)
    
                                        # tmp_TTS_add_pseudo = c(df_reads_filered$TTS, rep())
    if(nrow(df_reads_filered) == 0){
        return(NULL)
    }
                                        # stopifnot(nrow(df_reads_filered) > 0)
                                        # print(df_reads_filered)
    
    res_maxima = find_maxima(df_reads_filered$TTS, kde_bw = kde_bw, kde_min_gradient_length = kde_min_gradient_length)
    
    local_maxima = res_maxima$local_maxima
    global_maxima = res_maxima$global_maxima
    
                                        # global_count = calc_window_count(df_reads_labeled,
                                        #                                  window_center = global_maxima, 
                                        #                                  window_size = window_size)
    
                                        # local_counts = sapply(local_maxima,
                                        #                       function(x){
                                        #                         calc_window_count(df_reads_labeled,
                                        #                                           window_center = x, 
                                        #                                           window_size = window_size)
                                        #                       })
    
                                        # return(data.frame(
                                        #   pA = pA,
                                        #   strand = strand,
                                        #   local_maxima = paste0(local_maxima, collapse = ","),
                                        #   local_maxima_window_counts = paste0(local_counts, collapse = ","),
                                        #   global_maxima = global_maxima,
                                        #   global_maxima_window_count = global_count,
                                        #   total_count = nrow(df_reads_labeled)
                                        # ))
                                        # print(paste0("local maxima: ", paste0(local_maxima, collapse = ", ")))
                                        # print(paste0("global maxima: ", global_maxima))
                                        # 
    g1 = df_reads_labeled %>%
        ggplot(aes(x = TTS))+
        geom_histogram(bins = 100)+
        geom_vline(xintercept = local_maxima[local_maxima != global_maxima], color = "cyan")+
        geom_vline(xintercept = global_maxima, color = "hotpink")+
        labs(title = paste0("strand ",strand))+
                                        # scale_y_log10()+
                                        # facet_wrap(~ V14)+
        theme_bw()


    g2 = df_reads_labeled %>%
        ggplot(aes(x = TTS))+
        geom_histogram(bins = 100)+
        geom_vline(xintercept = local_maxima[local_maxima != global_maxima], color = "cyan")+
        geom_vline(xintercept = global_maxima, color = "hotpink")+
                                        # scale_y_log10()+
                                        # facet_wrap(~ V14)+
        theme_bw()

    if(length(global_maxima) > 0){
        g2 = g2+
            lims(x = c(global_maxima - 1000, global_maxima+1000))
    }
    return(g1 / g2)
}

                                        # 
                                        # tmp_reads = read_TLDR_reads_by_list_of_DHS(
                                        #   tmp_DHS$DHS,
                                        #   query_list_replicates = c("mES_WT_rep1", "mES_WT_rep2", "mES_WT_rep3",
                                        #                             "mES_WT_rep4", "mES_WT_rep5", "mES_WT_rep6")
                                        # )
                                        # 
                                        # find_peaks(df_reads = tmp_reads,
                                        #   pA = "+", 
                                        #   strand = tmp_DHS$strand_overlap_tss,
                                        #   kde_bw = 30,
                                        #   kde_min_gradient_length = 0,
                                        #   window_size = 30, # within a 30 bp window
                                        #   window_count_filter = 2) %>% 
                                        #   mutate(DHS = tmp_DHS$DHS,
                                        #          gene = tmp_DHS$gene_overlap_tss)
