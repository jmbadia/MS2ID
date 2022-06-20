#' @title Browse visually an Annot object
#'
#' @description
#' \code{MS2IDgui} is a Shiny app that visually browses the annotation results stored in an \linkS4class{Annot} object. The \linkS4class{Annot} object can be loaded as an argument or directly from the user interface.
#'
#' @details
#' The user interface has two panels. The \strong{right panel} lists the loaded
#' annotations as rows in a table; when one of them is clicked, the \strong{left
#' panel} responds by displaying the figures related to the selection:
#' \itemize{\item Tab QRY/REF: The figure compares the query (green) and the
#' reference (red) spectrum, with fragments in common (considering the
#' massErrorMSn argument) in a brighter colour; if present, precursor masses are
#' pointed out with a triangular symbol. Also, the bottom description shows the
#' inchikey of the reference compound and a link to a Pubchem's page listing
#' compounds with that inchikey.\item Tab Cons.: If the annotation has a
#' consensual query spectrum, this sub-panel displays a figure with the query
#' spectra that form it..}
#' The \strong{left panel} also has two more tabs with general features:
#' \itemize{\item Tab Home: Contains a control to load rds files (with
#' \linkS4class{Annot} object). \item Tab Info: Lists the arguments used in the
#' \code{\link{annotate}} function that resulted in the present
#' \linkS4class{Annot} object.}
#' Concerning the \strong{right panel}, it must be noted:
#' The \itemize{\item \code{mzML file} and \code{precursor M/z} controls filter
#' the annotations based in the query mzML file and the precursor mass of the
#' query spectrum, respectively. Also, the \code{redundant annotation} option
#' summarizes the annotations displaying only, for every query spectrum, the
#' best annotation per compound. \item Every row corresponds to an annotation in
#' the table, and different annotations of the same spectrum have identical row
#' colours. \item Column names that contain QRY, REF or CONS refer to query,
#' reference or consensus entities, respectively. e.g. \code{QRYprecursorMz}
#' contains the precursor M/z of the query spectra, \code{REFadduct} the adducts
#' of the reference spectra and \code{QRYrtime_CONS} the retention times of the
#' query spectra that conform every consensus spectrum. \item The column
#' \code{massNum} refers to the number of fragments of the query spectrum
#' (green), the reference spectrum (red) and the number of fragments in common
#' (black). Similarly, the column \code{collisionEnergy} shows the collision
#' energy of the query (green) and reference (red) spectra. \item
#' \code{REFinchikey} shows the inchikey of the reference compound and links to
#' a Pubchem's page with a list of compounds that shares that inchikey. \item
#' \code{ppmPrecMass} refers to the absolute difference value in ppm
#' between query and reference precursor masses.}
#'
#' @param annot \linkS4class{Annot} object with the results to be browsed. By
#'   default (no argument), the function opens the browser with no data and the
#'   \linkS4class{Annot} object can be loaded using the user interface.
#' @author Josep M. Badia \email{josepmaria.badia@@urv.cat}
#' @seealso \linkS4class{Annot} class and \code{\link{annotate}} function.
#' @example man/examples/loadMS2ID.R
#' @example man/examples/selectQuerySpectra.R
#' @examples
#'
#' ## ANNOTATION ---
#' library(MS2ID)
#' MS2IDlib <- MS2ID(MS2IDFolder)
#' annotResult <- annotate(QRYdata = queryFolder, MS2ID = MS2IDlib)
#'
#' ## Browse visually the results---
#' \dontrun{
#' MS2IDgui(annotResult)
#' }
#'
#' @import ggplot2
#' @import shiny
#' @importFrom shinyjs useShinyjs hidden hide toggle
#' @export
MS2IDgui <- function(annot){
    if(!missing(annot)){
        if(is(annot, "Annot")){
            .GlobalEnv$.jmb.rawAnnotSh <- annot
            nameAnnotObj <- as.list(match.call()[-1])$annot
            }else
            stop("'annot' argument is expected to be an Annot object")
    }
    #on.exit(rm(.jmb.rawAnnotSh, envir=.GlobalEnv))
    precDigits <- 4
    ui <- navbarPage(
        "MS2ID GUI",
        id = "tabs",
        theme = bslib::bs_theme(version = 4, bootswatch = "flatly"),
        #tabPanel("Home"),
        #insert here ui side of tabPanelIdentification.R, (in invisible2git)
        tabPanel(
            "Results",
            shinyjs::useShinyjs(),
            fluidPage(
                title = 'MS2ID GUI',
                fluidRow(
                    column(
                        4,
                        tabsetPanel(
                            type = "tabs",
                            tabPanel(
                                "",
                                class = "p-3 border border-top-0 rounded-bottom",
                                icon = icon("home", lib = "font-awesome"),
                                h4("Intro"), hr(),
                                htmlOutput("info"), br(),
                                h4("Settings"), hr(),
                                fileInput("lotFile", "Upload an Annot object",
                                          multiple = FALSE, accept = c(".rds")
                                          ),
                                ),
                            tabPanel(
                                "QRY/REF", br(),
                                h4(HTML(paste0(
                                    .customColor("QRY", "qry"), "/",
                                    .customColor("REF", "ref"), " MSn spectra")
                                    )),
                                shinyjs::hidden(
                                    div(
                                        id = "hiddenQRYREF",
                                        #content to show/hide
                                        HTML("<b><font color=\"#DF0054\">Pick a
                                             row to show the plot.</font></b>"))
                                ),
                                plotly::plotlyOutput('plotId'), br(),
                                htmlOutput("infoTabId")
                                ),
                            tabPanel(
                                "Cons", br(),
                                h4('QRY consensus formation'),
                                shinyjs::hidden(
                                    div(
                                        id = "hiddenCons",
                                        #content to show/hide
                                        HTML("<b><font color=\"#DF0054\">Pick a row to show the plot.</font></b>")),
                                    div(
                                        id = "hiddenIsCons",
                                        #content to show/hide
                                        HTML("<b><font color=\"#DF0054\">The selected spectrum is NOT a consensus spectrum.</font></b>"))
                                ),
                                plotly::plotlyOutput('plotCons'),
                                htmlOutput("infoTabCons")
                                ),
                            tabPanel(
                                "", br(),
                                h4("Annotation arguments"),

                                htmlOutput("infoTab"),
                                icon = icon("info", lib = "font-awesome")
                                )
                            ),
                        tags$hr(),
                        ),
                    column(
                        8,
                        htmlOutput("subheaderTable"),
                        tags$hr(),
                        fluidRow(
                            column(4, selectInput('ffile', 'mzML file', NULL)),
                            column(2, selectInput('UNKprec', 'precursor M/z ',
                                                  NULL, width=120)),
                            column(6, br(), br(),
                                   checkboxInput('reduntID',
                                                 'Redundant annotations ',
                                                 value = FALSE))
                        ),
                        DT::DTOutput('tbl'),
                        ),
                    ),
                )
            ),
        tabPanel("About",
                 HTML("
                 <h3>About</h3>
                 MS2IDgui is a Shiny app that visually browses the annotation results stored in an Annot object. Use the help available from the Help menu in R (help(MS2IDgui)) or visit <a href='https://jmbadia.github.io/MS2ID/articles/MS2ID.html'>the vignette</a> of the MS2ID package."))
        )


    server <- function(input, output, session) {
        observe({
            shinyjs::toggle(id = c("hiddenQRYREF"),
                            condition = !length(input$tbl_rows_selected))
            shinyjs::toggle(id = c("hiddenCons"),
                            condition = !length(input$tbl_rows_selected))
            if(!length(input$tbl_rows_selected)){
                shinyjs::hide(id = c("hiddenIsCons"))
                }
        })

        options(shiny.maxRequestSize=30*1024^2)
        #parent environment values
        colVisibles <- NA
        #insert here server side of tabPanelIdentification.R (in invisible2git)

        output$info <- renderText({.getFormText()})

        rawdata <- reactive({
            if(!is.null(input$lotFile)){
                lot <- readRDS(input$lotFile$datapath)
            }else if(exists(".jmb.rawAnnotSh", envir=.GlobalEnv))
                lot <- .GlobalEnv$.jmb.rawAnnotSh
            else
                req(FALSE)
            dt <- list(
                mtdt = .export2df(lot, summarizeHits = !input$reduntID),
                refSpctr = refSpectra(lot),
                qrySpctr = qrySpectra(lot),
                infoAnnot = infoAnnotation(lot)
                )
            dt$mtdt <- dt$mtdt %>%
                tidyr::replace_na(list(QRYdataOrigin="unknown",
                                       QRYprecursorMz="unknown"))
                dt$mtdt$massNum <- paste0(
                .customColor(dt$mtdt$QRYmassNum, "qry"), '/',
                dt$mtdt$cmnMasses,
                '/', .customColor(dt$mtdt$REFmassNum, "ref")
                )
            dt$mtdt$polarity <- paste0(
                .customColor(dt$mtdt$QRYpolarity, "qry"), '/',
                .customColor(dt$mtdt$REFpolarity, "ref")
                )
            dt$mtdt$collisionEnergy <- paste0(
                .customColor(dt$mtdt$QRYcollisionEnergy, "qry"), '/',
                .customColor(dt$mtdt$REFcollisionEnergy, "ref")
            )
            if(!"QRYacquisitionNum" %in% names(dt$mtdt)){
                dt$mtdt$QRYacquisitionNum <- NA
            }
            if("QRYacquisitionNum_CONS" %in% names(dt$mtdt) &
               all(is.na(dt$mtdt$QRYacquisitionNum))){
                dt$mtdt$QRYacquisitionNum <- dt$mtdt$QRYacquisitionNum_CONS
                }
            #round value to match with input select
            if("NLscore" %in% names(dt$mtdt))
                dt$mtdt$NLscore <- round(dt$mtdt$NLscore, 3)
            dt$mtdt$cosine <- round(dt$mtdt$cosine, 3)
            dt$mtdt <- .mergeCONS(dt$mtdt, fontsize = 1)
            #ORDER & SUBSET VISIBLE columns
            visiblVar <- c(
                "QRYrtime", "ppmPrecMass", INCRMETRIC, DECRMETRIC, "massNum",
                "propAdduct", "REFname", "REFformula", "REFadduct",
                "collisionEnergy", "REFpredicted", "REFID_db.spectra"
            )
            dt$mtdt <- dt$mtdt %>%
                select(-QRYmassNum, -cmnMasses, -REFmassNum, -QRYpolarity,
                       -REFpolarity, -QRYcollisionEnergy,
                       -REFcollisionEnergy) %>%
                #remove SOME empty columns
                select(propAdduct |
                           where(~ !(all(is.na(.)) | all(. == "")))) %>%
                select(sort(names(.))) %>%
                relocate(visiblVar[visiblVar %in% names(.)])
            #set columns visibility default
            colVisibles <<- names(dt$mtdt) %in% visiblVar
            return(dt)
        })

        metadataShowed <- reactive({
            req(input$UNKprec)
            rdMtdt <- rawdata()$mtdt
            if(is.numeric(rdMtdt$QRYprecursorMz))
                rdMtdt$QRYprecursorMz <- round(rdMtdt$QRYprecursorMz,
                                               precDigits)
            rdMtdt <- rdMtdt[
                isolate(rdMtdt$QRYdataOrigin == input$ffile) &
                    rdMtdt$QRYprecursorMz == input$UNKprec, ]
            #aesthetics modific.
            if("REFexactmass" %in% names(rdMtdt))
                rdMtdt$REFexactmass <- round(rdMtdt$REFexactmass, 3)
            return(rdMtdt)
        })

        getSelId <- reactive({
            rowSelected <- input$tbl_rows_selected
            if(!length(rowSelected)){
                req(FALSE)
            }
            mtdtShw <- isolate(metadataShowed()) %>%
                slice(rowSelected) %>%
                select(contains(c("idREFspect", "idQRYspect", "idREFcomp",
                                "QRYacquisitionNum","QRYrol")))
            return(mtdtShw)
            })

        output$tbl <- DT::renderDT(
            DT::datatable(
                metadataShowed(),
                caption = glue::glue("
                'file:' {isolate(input$ffile)}, QRY precursor m/z: \\
                {isolate(input$UNKprec)}
                    "),
                extensions = 'Buttons',
                options = list(
                    pageLength = 8,
                    info = FALSE,
                    lengthMenu = list(c(8,15,-1), c("8","15","All")),
                    columnDefs = list(list(visible=FALSE,
                                           targets=(which(!colVisibles)-1) )
                                      ),
                    stateLoadParams = DT::JS(
                        "function (settings, data) {return false;}"),
                    stateSave = TRUE,
                    scrollX = TRUE,
                    dom = 'Bfrt<"bottom"l>ip',
                    buttons = I('colvis')
                    ),
                selection = list(mode = 'single'),
                rownames= FALSE,
                escape = FALSE
            ) %>%
                DT::formatStyle(
                    'idQRYspect',
                    target = 'row',
                    backgroundColor = DT::styleEqual(
                        unique(metadataShowed()$idQRYspect),
                        head(
                            rep(c("#edf2fb","#ccdbfd"),
                                length(unique(metadataShowed()$idQRYspect))
                                ),
                            length(unique(metadataShowed()$idQRYspect))
                            )
                        )
                    )
            )

        output$subheaderTable <- renderText(
            if(exists(".jmb.rawAnnotSh", envir=.GlobalEnv)){
                paste(h3(paste0("Annotations in '", nameAnnotObj, "' variable"))
                      )
            }else{
                paste(h3("Please, load an Annot object"))
            }
        )

        output$infoTabId <- renderText({
            mtdtShw <- getSelId()
            vals2Show <- c(
                "cosine", "massNum", "propAdduct", "QRYprecursorMz", "REFname",
                "REFinchikey", "REFexactmass", "REFadduct", "REFprecursorMz",
                "REFpredicted", "REFID_db.comp")
            rd <- isolate(rawdata()$mtdt) %>%
                filter(idREFspect == mtdtShw$idREFspect,
                       idQRYspect == mtdtShw$idQRYspect,
                       idREFcomp == mtdtShw$idREFcomp) %>%
                select(contains(vals2Show))
            vals2Show[!vals2Show %in% names(rd)] <- "unknown"
            vals2Show[vals2Show %in% names(rd)] <- rd
            .getFormText("compare", as.list(vals2Show))
            #changed metadataShowed$REFinchikey[rowSelected] for metadata$REFinchikey
        })

        output$infoTab <- renderText({
            rdInfoAnnot <- isolate(rawdata()$infoAnnot)
            rdInfoAnnot$QRYdir <- rdInfoAnnot$... <- NULL
            rdInfoAnnot[rdInfoAnnot==""] <- "default"
            maxCharact <- 40
            toTrim <- rdInfoAnnot$QRYdata
            if(nchar(toTrim) > maxCharact + 3){
                rdInfoAnnot$QRYdata <- paste0(
                    "...",
                    substr(toTrim,
                           nchar(toTrim) - maxCharact + 1, nchar(toTrim)))
            }
            rdInfoAnnot$MS2ID <- basename(rdInfoAnnot$MS2ID)
            for(idArg in seq_along(rdInfoAnnot)){
                rdInfoAnnot$summ[[idArg]] <-
                    paste0("<b>", names(rdInfoAnnot)[idArg],
                           "</b>: ", rdInfoAnnot[[idArg]])
            }
            paste(rdInfoAnnot$summ, collapse = "<br> ")
        })

        output$infoTabCons <- renderText({
            mtdtShow <- getSelId()
            .getFormText("consens",
                         list(idQRY= mtdtShow$idQRYspect,
                              acqNum = mtdtShow$QRYacquisitionNum))
        })

        output$plotId <- plotly::renderPlotly({
            mtdtShow <- getSelId()
            rd <- rawdata()
            df1 <- .getSpectra2plot(rd$refSpctr, mtdtShow$idREFspect)
            df2 <- .getSpectra2plot(rd$qrySpctr, mtdtShow$idQRYspect)
            #only mz with hits
            mMz <- .matchMz(df2$x, df1$x,
                            isolate(rawdata()$infoAnnot$massErrMsn))
            hits1 <- !is.na(mMz)
            hits2 <- rep(FALSE, nrow(df2))
            hits2[mMz[hits1]] <- TRUE
            suppressWarnings(
                p <- ggplot() +
                    .get_gl(df1[!hits1,]) +
                    coord_cartesian(ylim = c(-100, 100)) +
                    .get_gl(df2[!hits2,], up = FALSE) +
                    .get_gl(df1[hits1,], up = TRUE, hit=TRUE) +
                    .get_gl(df2[hits2,], up = FALSE, hit=TRUE) +
                    labs(x = "m/z", y = "% intensity") +
                    theme_bw() + #black & white theme
                    .draw_precursor(df2, rd$mtdt, mtdtShow, nature = "qry") +
                    .draw_precursor(df1, rd$mtdt, mtdtShow, nature = "ref")
                )
            plotly::ggplotly(p, tooltip = c("text"))
        })

        output$plotCons <- plotly::renderPlotly({
            mtdtShow <- getSelId()
            shinyjs::toggle(id = c("hiddenIsCons"),
                            condition = (mtdtShow$QRYrol != 4))
            req(mtdtShow$QRYrol == 4)
            rd <- rawdata()
            df2 <- .getSpectra2plot(rd$qrySpctr, mtdtShow$idQRYspect)
            srcId <- unlist(rd$qrySpctr[rd$qrySpctr$id == mtdtShow$idQRYspect]$sourceSpect)
            vsd <- rd$qrySpctr[rd$qrySpctr$id %in% srcId]
            if("basePeakIntensity" %in% Spectra::spectraVariables(vsd)){
                vsdOrdr <- order(vsd$basePeakIntensity, decreasing = F)
            }else{
                vsdOrdr <- seq_along(nrow(vsd))
            }

            my_colors <- rep(.colCons, ceiling(length(vsd)/length(.colCons)))
            p <- ggplot() +
                geom_linerange(data = df2,
                               aes(x = x, ymax = y, ymin = 0,
                                   text=paste('</br>m/z: ',x,'</br>i: ',
                                              round(y,2))), colour="#00786C") +
                labs(x = "m/z", y = "% intensity") +
                theme_bw()#black & white theme

            for(src in vsdOrdr){
                df <- .getSpectra2plot(vsd[src])
                p <- p +
                    geom_linerange(data = df,
                                   aes(x = x, ymax = -y, ymin = 0,
                                       alpha=0.6,
                                       text=paste('</br>m/z: ', x,
                                                  '</br>i: ', round(y,2))
                                       ),
                                   colour = my_colors[src])
            }
            plotly::ggplotly(p, tooltip = c("text"))
        })

        observeEvent(input$lotFile, ignoreNULL = F, {
            freezeReactiveValue(input, "ffile")
            updateSelectInput(session, 'ffile',
                              choices = unique(rawdata()$mtdt$QRYdataOrigin)
            )
            if(exists(".jmb.rawAnnotSh", envir=.GlobalEnv)){
                req(input$lotFile)
                rm(.jmb.rawAnnotSh, envir=.GlobalEnv)
            }
            input$lotFile$name
            output$subheaderTable <- renderText(
                paste(h3(paste0("Annotations in '", input$lotFile$name,
                               "' file")))
            )
        })

        observeEvent(input$ffile, {
            freezeReactiveValue(input, "UNKprec")
            fileSel <- rawdata()$mtdt$QRYdataOrigin == input$ffile
            precSel <- unique(rawdata()$mtdt$QRYprecursorMz[fileSel])
            names <- suppressWarnings(as.numeric(precSel))
            if(!anyNA(precSel))
                names(precSel) <- round(names, precDigits)
            updateSelectInput(session, inputId = 'UNKprec',
                              choices = precSel)
            })

        observeEvent(input$UNKprec, {
            #save columns visibility before the inminnet table refresh
            if(!is.null(isolate(input$tbl_state$columns))) {
                colVisibles <<- isolate(
                    vapply(isolate(input$tbl_state$columns),
                           function(x) x$visible, FUN.VALUE = T)
                )}
        })
        }
    shinyApp(ui, server,  options = list(launch.browser = TRUE))
}
