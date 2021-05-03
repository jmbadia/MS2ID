#' @import ggplot2
#' @import shiny
MS2ID_gui <- function(){
    options(shiny.maxRequestSize=50*1024^2)

    ui <- fluidPage(
        title = 'MS2ID ID browser',
        fluidRow(
            column(4,
                   tabsetPanel(type = "tabs",
                               tabPanel("Options",hr(),
                                        fileInput("lotFile",
                                                  "Upload a lot file",
                                                  multiple = FALSE,
                                                  accept = c(".rds")),
                                        # Horizontal line ----
                                        tags$hr(),
                                        selectInput('ffile','Sample', NULL),
                                        selectInput('UNKprec',
                                                    'UNK precursor MZ ',
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

        output$info <- renderText({"Use the <font color=\"#3a86ff\">Options</font> tab to load a lot file (e.g. <a href=https://rovira-my.sharepoint.com/:u:/g/personal/39879942-n_epp_urv_cat/ERSapzHJPIBNgTxxYYKqqZsBMoWKj3ag99l-TqZOXDK43Q?e=PchTkL target='_blank'>this sample</a>) and configure. Once the Identifications table is loaded, click a row to see the QRY/REF MS2 spectra on the <font color=\"#3a86ff\">Plot</font> tab"})

        rawdata <- reactive({
            req(input$lotFile)
            lot <- readRDS(input$lotFile$datapath)
            dt <- list(mtdt = .export2df(lot), refSpctr = refSpectra(lot),
                       qrySpctr = qrySpectra(lot))
            #lotMetdt <- .export2df(lot)
            dt$mtdt$massNum <- paste0('<font color=\"#00786C\">',
                                       dt$mtdt$QRYmassNum,'</font>/',
                                       dt$mtdt$cmnMasses,'/<font color=\"#d64c1d\">',
                                       dt$mtdt$REFmassNum,'</font>')
            #obtain polarity
            dt$mtdt$polarity <- paste0('<font color=\"#00786C\">',
                                        dt$mtdt$QRYpolarity,
                                        '</font>/<font color=\"#d64c1d\">',
                                        dt$mtdt$REFpolarity,
                                        '</font>')
            dt$mtdt$CE <- paste0('<font color=\"#00786C\">', dt$mtdt$QRYcollisionEnergy,
                                  '</font>/<font color=\"#d64c1d\">',
                                  dt$mtdt$REFcollisionEnergy, '</font>')

            if(all(is.na(dt$mtdt$QRYacquisitionNum)))
                dt$mtdt$QRYacquisitionNum <- dt$mtdt$QRYacquisitionNum_CONS

            #round value to match with input select
            dt$mtdt$cosine <- round(dt$mtdt$cosine, 3)

            #remove redundant info
            dt$mtdt$QRYacquisitionNum_CONS <- dt$mtdt$QRYmassNum <-
                dt$mtdt$cmnMasses <- dt$mtdt$REFmassNum <- NULL

            dt$mtdt <- dplyr::relocate(dt$mtdt, CE, .after=QRYrtime)
            dt$mtdt <- dplyr::relocate(dt$mtdt, massNum, .after=REFformula)
            dt$mtdt <- dplyr::relocate(dt$mtdt, QRYacquisitionNum, .before=QRYrtime)
            #ORDER & SUBSET columns
            visiblVar <- c("QRYrtime", "QRYacquisitionNum",
                         "REFMmi","propAdduct","REFadduct",
                         "REFprecursorMz", INCRMETRIC, DECRMETRIC, "REFname","REFformula",
                         "massNum", "CE","QRYprecInt",
                         "REFnature","REFinstrument",
                         "REFinstrumentType","idQRYspect",
                         "idREFspect","idREFcomp","REFID_db.comp", "REFID_db.spect")

            #First time, set columns visibility default
            if(is.na(colVisibles))
                colVisibles <<- names(dt$mtdt) %in% visiblVar
            return(dt)
        })

        metadataShowed <- reactive({
            req(input$UNKprec)
            #rdmetadata <- rawdata()$metadata
            rdmetadata <- rawdata()$mtdt
            #round value to match with input select
            rdmetadata$QRYprecursorMz <- round(rdmetadata$QRYprecursorMz, 4)
            rdmetadata <- rdmetadata[isolate(rdmetadata$QRYdataOrigin==input$ffile) &
                                         rdmetadata$QRYprecursorMz==input$UNKprec,]
            #select QRYacq according if they are consensus or not
            #browser()

            if(!input$reduntID){ #if NOT checked, keep only the best score of every scan-REFmetabolite
                rdmetadata <- rdmetadata %>%
                    group_by(idQRYspect, idREFcomp) %>%
                    top_n(1, cosine) %>%
                    distinct(idQRYspect, idREFcomp, .keep_all = T) %>%
                    ungroup() %>% arrange(QRYprecursorMz, idQRYspect) %>%
                    as.data.frame()
            }
            #round leftovers
            rdmetadata$REFMmi <- round(rdmetadata$REFMmi, 4)
            rdmetadata$QRYrtime <- round(rdmetadata$QRYrtime, 2)
            rdmetadata$QRYacquisitionNum <- paste0("<font size=1>", rdmetadata$QRYacquisitionNum,"</font>")
            return(rdmetadata)
        })

        output$tbl <- DT::renderDT(
            DT::datatable(metadataShowed(),
                      caption = paste('file:',isolate(input$ffile),
                                      ', QRY precursor m/z:',
                                      isolate(as.numeric(input$UNKprec))),
                      colnames = c('QRYprecMz' = 'QRYprecursorMz',
                                   'REFprecMz' = 'REFprecursorMz',
                                   'QRYacqNum' = 'QRYacquisitionNum'),
                      extensions = 'Buttons',
                      options = list(
                          pageLength = 8, info = FALSE,
                          lengthMenu = list(c(8,15,-1), c("8","15","All")),
                          columnDefs = list(list(visible=FALSE,
                                                 targets=(which(!colVisibles)-1) )),
                          stateLoadParams = DT::JS("function (settings, data) {return false;}"),
                          stateSave=TRUE,
                          scrollX = TRUE,
                          dom = 'Bfrt<"bottom"l>ip',
                          buttons = I('colvis')),
                      selection = list(mode = 'single'),
                      rownames= FALSE,
                      escape = FALSE
            )%>% DT::formatStyle(
                'QRYacqNum',
                target = 'row',
                #QRYacquisitionNum_CONS <- idQRYspect (3)
                backgroundColor = DT::styleEqual(unique(metadataShowed()$idQRYspect),
                                             head(rep(c("#f8f8f8","#eae9e9"),
                                                      length(unique(metadataShowed()$idQRYspect))),
                                                  length(unique(metadataShowed()$idQRYspect)))
                ))
            )

        output$subheaderTable <- renderText(
            paste(" <link rel='stylesheet' href='https://fonts.googleapis.com/css?family=Montserrat'><span style='font-size: 24pt;font-weight:bold;font-family: Montserrat'>Identifications </span>")
        )

        #To show value live
        #output$esborra <-  renderText({ })

        output$REFmetadata <- renderText({
            rowSelected <- input$tbl_rows_selected
            if(!length(rowSelected)){
                return()
            }
            metadataShowed <- isolate(metadataShowed())
            rd <- isolate(rawdata()$mtdt)
            metadata <- rd[rd$idREFspect == metadataShowed[rowSelected, "idREFspect"] &
                               rd$idQRYspect==metadataShowed[rowSelected, "idQRYspect"] &
                               rd$idREFcomp==metadataShowed[rowSelected, "idREFcomp"],]
            paste0("<div style='text-align:center;'><font color=\"#333333\"><br>","cosine= ",
                   metadata$cosine  ,", massNum= ", metadata$massNum, "<br></font> </div>",
                   "<font color=\"#00786C\"><font size=4><b>QRY</b></font><br>",
                   "<b>propAdduct</b>= ", metadata$propAdduct,", <b>precursor</b>= ",
                   metadata$QRYprecursorMz,"</font><br>",
                   "<font color=\"#d64c1d\"><font size=4><b>REF</b><br></font>",
                   metadata$REFname,
                   '<a href="https://www.ncbi.nlm.nih.gov/pccompound?term=%22',
                   metadataShowed$REFinchikey[rowSelected], '%22[InChIKey]" title="link to Pubchem using REFinchiKey as search word" target="_blank"> (PubChem)</a>',
                   "<br><b>Mmi</b>= ", metadata$REFMmi,", <b>adduct</b>= ", metadata$REFadduct,", <b>precursor</b>= ", metadata$REFprecursorMz,
                   ", <br><b>nature=</b> ", metadata$REFnature,", <b>DB=</b> ", metadata$REFID_db.comp,"</font>")
        })

        output$plot <- plotly::renderPlotly({
            rowSelected <- input$tbl_rows_selected
            if(!length(rowSelected)){
                return()
            }
            metadataShowed <- isolate(metadataShowed())
            rd <- rawdata()
            refSpectra <- rd$refSpctr[rd$refSpctr$id == metadataShowed[rowSelected, "idREFspect"]]
            refSp_mz <- unlist(Spectra::mz(refSpectra))
            refSp_i <- unlist(Spectra::intensity(refSpectra))
            refSp_i <- 100*refSp_i/max(refSp_i)

            unkSpectra <- rd$qrySpctr[rd$qrySpctr$id == metadataShowed[rowSelected, "idQRYspect"]]
            qrySp_mz <- unlist(Spectra::mz(unkSpectra))
            qrySp_i <- unlist(Spectra::intensity(unkSpectra))
            qrySp_i <- 100*qrySp_i/max(qrySp_i)

            df1<-data.frame(x = refSp_mz, y = refSp_i)
            df2<-data.frame(x = qrySp_mz, y = qrySp_i)

            p <- ggplot() +
                geom_linerange(data = df1,
                               aes(x = x, ymax=-y, ymin=0,
                                   text=paste('</br>m/z: ',x,'</br>i: ',
                                              round(y,2))), colour="#d64c1d") +
                coord_cartesian(ylim = c(-100, 100))
            p <- p +
                geom_linerange(data = df2,
                               aes(x = x, ymax = y, ymin = 0,
                                   text=paste('</br>m/z: ',x,'</br>i: ',
                                              round(y,2))), colour="#00786C")
            p <- p +
                labs(title = "QRY/REF MS2 spectra", x="m/z", y="% intensity")
            p <- p + theme_bw()#black & white theme
            #add precursor representation
            ##obtain precursors
            UNKprec_mz <- unique(rd$mtdt[rd$mtdt$idQRYspect ==
                                             metadataShowed[rowSelected, "idQRYspect"],"QRYprecursorMz"])
            if(!is.na(UNKprec_mz)){
                UNKprec_int <- 100 #default
                near_prec <- abs(qrySp_mz - UNKprec_mz) < 0.01
                if(any(near_prec))
                    UNKprec_int <- 15 + max(qrySp_i[near_prec])
                p <- p + geom_point(aes(y = UNKprec_int, x = UNKprec_mz),
                                    shape = 25, colour="#00786C")
            }
            REFprec_mz <- unique(rd$mtdt[rd$mtdt$idREFspect == metadataShowed[rowSelected, "idREFspect"],"REFprecursorMz"])
            if(!is.na(REFprec_mz)){
                REFprec_int <- 100 #default
                near_prec <- abs(refSp_mz - REFprec_mz) < 0.01
                if(any(near_prec))
                    REFprec_int <- 15 + max(refSp_i[near_prec])
                p <- p + geom_point(aes(y = -REFprec_int, x = REFprec_mz),
                                    shape = 24, colour="#d64c1d")
            }
            plotly::ggplotly(p, tooltip = c("text"))
        })

        observeEvent(input$lotFile, {
            updateSelectInput(session,'ffile',
                              choices=unique(rawdata()$mtdt$QRYdataOrigin)
            )
        })

        observeEvent(input$ffile, {
            updateSelectInput(session,'UNKprec',
                              choices=round(unique(rawdata()$mtdt$QRYprecursorMz[rawdata()$mtdt$QRYdataOrigin==input$ffile]),4))
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
