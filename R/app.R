#' @import ggplot2
MS2ID_gui <- function(){
    options(shiny.maxRequestSize=50*1024^2)

    ui <- fluidPage(
        title = 'MS2ID ID browser',
        fluidRow(
            column(4,
                   tabsetPanel(type = "tabs",
                               tabPanel("Options",hr(),
                                        fileInput("lotFile", "Upload a lot file",
                                                  multiple = FALSE,
                                                  accept = c(".rds")),
                                        # Horizontal line ----
                                        tags$hr(),
                                        selectInput('ffile','Sample', NULL),
                                        selectInput('UNKprec','UNK precursor MZ ',
                                                    NULL),
                                        checkboxInput('reduntID',
                                                      'Show redundant identifications ',
                                                      value = FALSE),
                               ),
                               tabPanel("Plot", hr(),
                                        plotly::plotlyOutput('plot'), br(),
                                        htmlOutput("REFmetadata"),
                               ),
                               tabPanel("INFO", hr(), htmlOutput("info"))
                   ),
                   tags$hr(),
            ),
            column(8,
                   htmlOutput("subheaderTable"),
                   tags$hr(),
                   DT::DTOutput('tbl'),
            ),
        ),  )


    server <- function(input, output, session) {
        #parent environment values
        colVisibles <- NA

        output$info <- renderText({"Use the <font color=\"#3a86ff\">Options</font> tab to load a lot file (e.g. <a href=https://rovira-my.sharepoint.com/:u:/g/personal/39879942-n_epp_urv_cat/ERSapzHJPIBNgTxxYYKqqZsBMoWKj3ag99l-TqZOXDK43Q?e=PchTkL target='_blank'>this sample</a>) and configure. Once the Identifications table is loaded, click a row to see the UNK/REF MS2 spectra on the <font color=\"#3a86ff\">Plot</font> tab"})

        rawdata <- reactive({
            req(input$lotFile)
            lot <- readRDS(input$lotFile$datapath)

            #if we dont have REF mass number, we calculate it. REMOVE THIS, mtrust we always will have the REF mass number.
            if(!"REFmassNum" %in% names(lot$metadata)){
                #obtain mass number on REF
                REFmassNum <- vapply(lot$REFfragments$spectra, ncol, FUN.VALUE = 2)
                REFmassNum <- REFmassNum[match(lot$metadata$REFidSpectra, lot$REFfragments$ID_spectra)]
                lot$metadata$REFmassNum <- REFmassNum
            }
            lot$metadata$massNum <- paste0('<font color=\"#00786C\">',lot$metadata$UNKmassNum,'</font>/',lot$metadata$MATCHmassNum,'/<font color=\"#d64c1d\">',lot$metadata$REFmassNum,'</font>')
            #obtain polarity
            lot$metadata$polarity <- paste0('<font color=\"#00786C\">', lot$metadata$UNKpolarity, '</font>/<font color=\"#d64c1d\">', lot$metadata$REFpolarity,'</font>')
            #obtain CE
            lot$metadata$CE <- paste0('<font color=\"#00786C\">', lot$metadata$UNK_CE, '</font>/<font color=\"#d64c1d\">', lot$metadata$REF_CE,'</font>')
            #SUBSET & reorder cols
            lot$metadata <- lot$metadata[,c("UNKprecMZ", "UNKrt", "cossim","massNum", "assmdUNKAdduct", "REFname","REFformula", "REFMmi", "polarity","CE","REFidDB" ,"UNKacqNum","UNKprecCharge","UNKprecInt","REFnature","REFinstrument","REFionSource","REFcasNum","REFprecMZ","REFprecAdduct","REFinchikey","UNKmsLevel","UNKidSpectra","REFidSpectra","REFidMetabolite","file")]

            #First time, set columns visibility default
            if(is.na(colVisibles)) colVisibles <<- !(names(lot$metadata) %in% c("file","UNKacqNum","UNKprecMZ","REFprecAdduct","REFprecMZ","REFinchikey","REFcasNum","UNKprecCharge","UNKprecInt","REFnature","REFinstrument","REFionSource","UNKmsLevel","UNKidSpectra","REFidSpectra","REFidMetabolite","REFidDB"))
            return(lot)
        })

        metadataShowed <- reactive({
            req(input$UNKprec)
            rdmetadata <- rawdata()$metadata
            #round value to match with input select
            rdmetadata$UNKprecMZ <- round(rdmetadata$UNKprecMZ,4)
            rdmetadata <- rdmetadata[isolate(rdmetadata$file==input$ffile) & rdmetadata$UNKprecMZ==input$UNKprec,]

            if(!input$reduntID){ #if NOT checked, keep only the best score of every scan-REFmetabolite
                rdmetadata <- rdmetadata %>%
                    dplyr::group_by(UNKidSpectra, REFidMetabolite) %>%
                    dplyr::top_n(1, cossim) %>%
                    dplyr::distinct(UNKidSpectra, REFidMetabolite, .keep_all = T) %>%
                    dplyr::ungroup() %>%
                    dplyr::arrange(UNKprecMZ, UNKidSpectra) %>%
                    as.data.frame()
            }
            #round leftovers
            rdmetadata$REFMmi <- round(rdmetadata$REFMmi,4)
            rdmetadata$UNKrt <- round(rdmetadata$UNKrt,2)
            return(rdmetadata)
        })

        output$tbl <- DT::renderDT(
            DT::datatable(metadataShowed(),
                      caption = paste('file:',isolate(input$ffile),', UNK precursor m/z:', isolate(as.numeric(input$UNKprec))),
                      colnames = c('UNKaddct?' = 'assmdUNKAdduct'),
                      extensions = 'Buttons',
                      options = list(
                          pageLength = 8, info = FALSE,
                          lengthMenu = list(c(8,15,-1), c("8","15","All")),
                          columnDefs = list(list(visible=FALSE, targets=(which(!colVisibles)-1) )),
                          stateLoadParams = DT::JS("function (settings, data) {return false;}"),
                          stateSave=TRUE,
                          scrollX = TRUE,
                          dom = 'Bfrt<"bottom"l>ip',
                          buttons = I('colvis')),
                      selection = list(mode = 'single'),
                      rownames= FALSE,
                      escape = FALSE
            )%>% DT::formatStyle(
                'UNKacqNum',
                target = 'row',
                backgroundColor = DT::styleEqual(unique(metadataShowed()$UNKacqNum),
                                             head(rep(c("#f8f8f8","#eae9e9"),length(unique(metadataShowed()$UNKacqNum))),length(unique(metadataShowed()$UNKacqNum)))
                )))

        output$subheaderTable <- renderText(
            paste(" <link rel='stylesheet' href='https://fonts.googleapis.com/css?family=Montserrat'><span style='font-size: 24pt;font-weight:bold;font-family: Montserrat'>Identifications </span>")
        )

        #To show value live
        #output$esborra <-  renderText({ })

        output$REFmetadata <-  renderText({
            rowSelected <- input$tbl_rows_selected
            if(!length(rowSelected)){
                return()
            }
            metadataShowed <- isolate(metadataShowed())
            rd <- isolate(rawdata())
            metadata <- rd$metadata[rd$metadata$REFidSpectra==metadataShowed[rowSelected, "REFidSpectra"] & rd$metadata$UNKidSpectra==metadataShowed[rowSelected, "UNKidSpectra"] & rd$metadata$REFidMetabolite==metadataShowed[rowSelected, "REFidMetabolite"],]
            paste0("<div style='text-align:center;'><font color=\"#333333\"><br>","cossim= ", metadata$cossim  ,", massNum= ", metadata$massNum, "<br></font> </div>",
                   "<font color=\"#00786C\"><font size=4><b>UNK</b></font><br>", "<b>UNKaddct?</b>= ", metadata$assmdUNKAdduct,", <b>precMZ</b>= ", metadata$UNKprecMZ,"</font><br>",
                   "<font color=\"#d64c1d\"><font size=4><b>REF</b><br></font>", metadata$REFname,
                   '<a href="https://www.ncbi.nlm.nih.gov/pccompound?term=%22',metadataShowed$REFinchikey[rowSelected], '%22[InChIKey]" title="link to Pubchem using REFinchiKey as search word" target="_blank"> (PubChem)</a>',
                   "<br><b>Mmi</b>= ", metadata$REFMmi,", <b>adduct</b>= ", metadata$REFprecAdduct,", <b>precMZ</b>= ", metadata$REFprecMZ,
                   ", <br><b>nature=</b> ", metadata$REFnature,", <b>DB=</b> ", metadata$REFidDB,"</font>")
        })

        output$plot <- plotly::renderPlotly({
            rowSelected <- input$tbl_rows_selected
            if(!length(rowSelected)){
                return()
            }
            metadataShowed <- isolate(metadataShowed())
            rd <- rawdata()
            refSpectra <- rd$REFfragments$spectra[[which(rd$REFfragments$ID_spectra == metadataShowed[rowSelected, "REFidSpectra"])]]
            refSpectra["intensity",] <- 100*refSpectra["intensity",]/max(refSpectra["intensity",])
            unkSpectra <- rd$UNKfragments$spectra[[which(rd$UNKfragments$idUNKSpectra == metadataShowed[rowSelected, "UNKidSpectra"])]]
            unkSpectra["intensity",] <- 100*unkSpectra["intensity",]/max(unkSpectra["intensity",])
            df1<-data.frame(x = refSpectra[1,], y = refSpectra[2,])
            df2<-data.frame(x = unkSpectra[1,], y = unkSpectra[2,])

            p <- ggplot() +
                geom_linerange(data = df1, aes(x = x, ymax=-y, ymin=0, text=paste('</br>m/z: ',x,'</br>i: ', round(y,2))), colour="#d64c1d") +
                coord_cartesian(ylim = c(-100, 100))
            p <- p + geom_linerange(data = df2, aes(x = x, ymax = y, ymin = 0, text=paste('</br>m/z: ',x,'</br>i: ',round(y,2))), colour="#00786C")
            p <- p + labs(title = "UNK/REF MS2 spectra", x="m/z", y="% intensity")
            p <- p + theme_bw()#black & white theme
            #add precursor representation
            ##obtain precursors
            UNKprec_mz <- unique(rd$metadata[rd$metadata$UNKidSpectra == metadataShowed[rowSelected, "UNKidSpectra"],"UNKprecMZ"])
            if(!is.na(UNKprec_mz)){
                UNKprec_int <- 100 #default
                near_prec <- abs(unkSpectra[1,]-UNKprec_mz)<0.01
                if(any(near_prec))
                    UNKprec_int <- 15 + max(unkSpectra[2, near_prec])
                p <- p + geom_point(aes(y = UNKprec_int, x = UNKprec_mz),
                                    shape = 25, colour="#00786C")
            }
            REFprec_mz <- unique(rd$metadata[rd$metadata$REFidSpectra == metadataShowed[rowSelected, "REFidSpectra"],"REFprecMZ"])
            if(!is.na(REFprec_mz)){
                REFprec_int <- 100 #default
                near_prec <- abs(refSpectra[1,]-REFprec_mz)<0.01
                if(any(near_prec)) REFprec_int <- 15 + max(refSpectra[2, near_prec])
                p <- p + geom_point(aes(y = -REFprec_int, x = REFprec_mz), shape = 24, colour="#d64c1d")
            }
            plotly::ggplotly(p, tooltip = c("text"))
        })

        observeEvent(input$lotFile, {
            updateSelectInput(session,'ffile',
                              choices=unique(rawdata()$metadata$file)
            )
        })

        observeEvent(input$ffile, {
            updateSelectInput(session,'UNKprec',
                              choices=round(unique(rawdata()$metadata$UNKprecMZ[rawdata()$metadata$file==input$ffile]),4))
        })

        observeEvent(input$UNKprec, {
            #save columns visibility before the inminnet table refresh
            if(!is.null(isolate(input$tbl_state$columns))) {
                colVisibles <<- isolate(vapply(isolate(input$tbl_state$columns), function(x) x$visible, FUN.VALUE = T))
            }
        })
    }
    shinyApp(ui,server)
}
