.getFormText <- function(type = "info", varList){
    switch (type,
        info = glue::glue("
        Use the <i>Settings</i> section to load a rds file with an Annot \\
        object (e.g. AnnotSample.rds hosted in \\
        <a href=https://doi.org/10.5281/zenodo.5561324 target='_blank'> \\
        Zenodo</a>). Once the \\
        annotations table (right panel) is loaded, select a row (query
        spectrum) and check its QRY vs REF spectra plot (<i>\\
        QRY/REF</i> tab) or its consensus formation (<i>Cons.</i> tab).
                          "),
        consens = glue::glue("
        <br><font color=\"#00786C\"><b>This query consensus spectrum</b> \\
        ({varList[[1]]})</font> has been obtained from query spectra with the \\
        acquisition numbers {varList[[2]]}.
                          "),
        compare = glue::glue("
        <div style=\"text-align:center;\"><font color=\"#333333\"><br>cosine= \\
        {varList[[1]]}, massNum= {varList[[2]]}<br></font> </div><font \\
        color=\"#00786C\"><font size=4><b>QueRY spectra</b></font>\\
        <br><b>propAdduct</b>={varList[[3]]}, <b>precursor</b>= {varList[[4]]}\\
        </font><br><font color=\"#d64c1d\"><font size=4><b>\\
        REFerence spectra</b><br></font>{varList[[5]]}<a \\
        href=\"https://www.ncbi.nlm.nih.gov/pccompound?term=%22{varList[[6]]}\\
        %22[InChIKey]\" title=\"link to Pubchem using REFinchiKey as search \\
        word\" target=\"_blank\"> (PubChem)</a><br><b>Mmi</b>= {varList[[7]]},\\
        <b>adduct</b>= {varList[[8]]}, <b>precursor</b>= {varList[[9]]}, <br>\\
        <b>nature=</b> {varList[[10]]}, <b>DB=</b>{varList[[11]]}</font>
                          ")
    )
}

#obtain the spectrum that matches the id in a df form in order 2 plot it
.getSpectra2plot <- function(rawSpectr, id){
    if(!missing(id)){
        rawSpectr <- rawSpectr[rawSpectr$id == id]
    }
    Sp_mz <- unlist(Spectra::mz(rawSpectr))
    Sp_i <- unlist(Spectra::intensity(rawSpectr))
    Sp_i <- 100*Sp_i/max(Sp_i)
    return(data.frame(x = Sp_mz, y = Sp_i))
}

.get_gl <- function(df, up = TRUE, hit = FALSE){
    if(nrow(df) == 0) return(NULL)
    a <- ifelse(up, -1, 1)
    colour <- ifelse(hit,
                     ifelse(up, .MS2IDcols["ref"], .MS2IDcols["qry"]),
                     ifelse(up,
                            .MS2IDcols["refLight"], .MS2IDcols["qryLight"]))
    geom_linerange(
        data = df,
        aes(
            x = x,
            ymax = a*y,
            ymin = 0,
            text = paste('</br>m/z: ', x ,'</br>i: ', round(y, 2))
        ),
        colour = colour
    )
}

.draw_precursor <- function(df, mtdt, mtdtShw, nature = "qry"){
    if(nrow(df) == 0) return(NULL)
    if(nature == "qry" & "QRYprecursorMz" %in% names(mtdt)){
        a = 1
        shape = 25
        colour = "#00786C"
        prec_mz <- mtdt %>%
            filter(idQRYspect == mtdtShw$idQRYspect) %>%
            select(QRYprecursorMz)
    }else if(nature == "ref" & "REFprecursorMz" %in% names(mtdt)){
        a = -1
        shape = 24
        colour = "#d64c1d"
        prec_mz <- mtdt %>%
            filter(idREFspect == mtdtShw$idREFspect) %>%
            select(REFprecursorMz)
    }else
        return()
    prec_mz <- as.numeric(unique(unlist(prec_mz)))

    if(!is.na(prec_mz) & prec_mz != "unknown"){
        prec_int <- 100 #default
        near_prec <- abs(df$x - prec_mz) < 0.01
        if(any(near_prec)){
            prec_int <- 15 + max(df$y[near_prec])
        }
        geom_point(
            aes(y = a * prec_int, x = prec_mz),
            shape = shape, colour = colour
        )
    }else{
        geom_point()
    }
}

.customColor <- function(txt, type){
    initFont <- switch(type, ref = .MS2IDcols["ref"], qry = .MS2IDcols["qry"])
    paste0('<font color=\"', initFont, '\">', txt, '</font>')
}
