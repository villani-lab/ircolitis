"use strict";

var g_bins = null

var pals = {
  "alphabet": [
    "#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919",
    "#005C31", "#2BCE48", "#FFCC99", "#808080", "#94FFB5", "#8F7C00",
    "#9DCC00", "#C20088", "#003380", "#FFA405", "#FFA8BB", "#426600",
    "#FF0010", "#5EF1F2", "#00998F", "#E0FF66", "#740AFF", "#990000",
    "#FFFF80", "#FFE100", "#FF5005",
    "#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919",
    "#005C31", "#2BCE48", "#FFCC99", "#808080", "#94FFB5", "#8F7C00",
    "#9DCC00", "#C20088", "#003380", "#FFA405", "#FFA8BB", "#426600",
    "#FF0010", "#5EF1F2", "#00998F", "#E0FF66", "#740AFF", "#990000",
    "#FFFF80", "#FFE100", "#FF5005"
  ],
  "okabe": [
    "#666666", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7"
  ],
  "mpn65": [
    '#ff0029', '#377eb8', '#66a61e', '#984ea3', '#00d2d5', '#ff7f00', '#af8d00',
    '#7f80cd', '#b3e900', '#c42e60', '#a65628', '#f781bf', '#8dd3c7', '#bebada',
    '#fb8072', '#80b1d3', '#fdb462', '#fccde5', '#bc80bd', '#ffed6f', '#c4eaff',
    '#cf8c00', '#1b9e77', '#d95f02', '#e7298a', '#e6ab02', '#a6761d', '#0097ff',
    '#00d067', '#000000', '#252525', '#525252', '#737373', '#969696', '#bdbdbd',
    '#f43600', '#4ba93b', '#5779bb', '#927acc', '#97ee3f', '#bf3947', '#9f5b00',
    '#f48758', '#8caed6', '#f2b94f', '#eff26e', '#e43872', '#d9b100', '#9d7a00',
    '#698cff', '#d9d9d9', '#00d27e', '#d06800', '#009f82', '#c49200', '#cbe8ff',
    '#fecddf', '#c27eb6', '#8cd2ce', '#c4b8d9', '#f883b0', '#a49100', '#f48800',
    '#27d0df', '#a04a9b'
  ],
  "batlow": [
    "#001959", "#021C5A", "#041F59", "#04215B", "#06245B", "#07265B",
    "#09295C", "#0A2C5C", "#0A2E5C", "#0C315D", "#0C335D", "#0D365D",
    "#0F395F", "#0F3C5F", "#103F60", "#114160", "#124460", "#134660",
    "#154A61", "#164C60", "#184F60", "#1A5260", "#1C5460", "#1E5761",
    "#20595F", "#235C5F", "#255E5E", "#29605E", "#2C635C", "#2F655B",
    "#336759", "#366858", "#3A6A56", "#3E6C54", "#416D53", "#456F50",
    "#49704F", "#4D724C", "#50734B", "#547548", "#597645", "#5C7744",
    "#607841", "#647940", "#687B3D", "#6C7B3C", "#707D39", "#757E37",
    "#797F35", "#7D8133", "#818132", "#868330", "#8A842F", "#90862E",
    "#95872D", "#9A872D", "#9F892D", "#A48A2D", "#AA8C2F", "#AF8D30",
    "#B58D33", "#BB8F36", "#C08F38", "#C5903C", "#CA913F", "#CF9243",
    "#D49347", "#D9944C", "#DE9651", "#E29755", "#E6985A", "#E9995F",
    "#ED9B65", "#F19D6B", "#F39E70", "#F5A076", "#F7A17B", "#F9A380",
    "#F9A486", "#FBA68B", "#FBA892", "#FDA997", "#FCAC9C", "#FDADA1",
    "#FDAFA6", "#FDB0AC", "#FDB3B2", "#FDB5B7", "#FDB6BC", "#FDB9C2",
    "#FDBAC7", "#FDBCCC", "#FDBDD2", "#FCC0D8", "#FBC2DD", "#FBC3E3",
    "#FAC6E8", "#FAC7EE", "#FBCAF4", "#F9CCF9"
  ].reverse(),
}

// Store one gene in a list of objects.
var g_mydata = []
var g_font_size = 12

var g_healthField = "case"
var g_donorField = "donor"

var g_count = 0
var g_loadGeneAndColor = 0
var g_colorByDefault = 0

var db = null
var g_coords = null
var g_exprArr, g_decArr, g_geneSym, g_geneDesc, g_binInfo

var g_blob = null

var g_size = {
  gene_bars: {
    width: 160,
    height: (n_groups, n_subgroups) => Math.max(20 * n_groups * (0.5 * n_subgroups), 300)
  },
  gene_boxplot: {
    width: 200,
    height: (n_groups, n_subgroups) => Math.max(20 * n_groups * (0.5 * n_subgroups), 300),
    margin: {
      top: 15, right: 1, bottom: 40, left: 5, xoffset: 160
    }
  },
  hex_radius: 3
}

var mybrowser = function() {
  // var db = null; // the cbData object from cbData.js. Loads coords,

  var gRecentGenes = [];
  // var coords = null;

  // var vlSpec = vega_geneHeatmap;

  var recent_genes = []

  var state = {
    ds: "n3_2",
    gene: "none",
    gene_groupby: "none",
    gene_loaded: false,
    meta_colorby: "cluster",
    query_geneinfo: false
  }

  $("body").html(
    `
  <div class="navbar navbar-expand-md navbar-light bg-light py-4">
    <div class="container">
      <a class="navbar-brand" href="/">Villani Lab</a>
      <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarResponsive" aria-controls="navbarResponsive" aria-expanded="false" aria-label="Toggle navigation">
        <span class="navbar-toggler-icon"></span>
      </button>
      <div class="collapse navbar-collapse sidebar-nav" id="navbarResponsive">
        <ul class="navbar-nav ml-auto mr-auto">
          <li class="nav-item"><a class=
          "nav-link mr-3 mt-2 mb-1" href="/ircolitis" title=
          "irColitis">irColitis</a></li>
          <li class="nav-item active"><a class=
          "nav-link mr-3 mt-2 mb-1" href="/ircolitis/app/?ds=all_cells_data&gene=RORA&groupby=none"
          title="Cell Clusters">Cell Clusters</a></li>
          <li class="nav-item"><a id="link-gene-contrasts" class=
          "nav-link mr-3 mt-2 mb-1" href="/ircolitis/gene-contrasts/"
          title="Gene Contrasts">Gene Contrasts</a></li>
        </ul>
      </div>
    </div>
  </div>

    <div class="container" id="mycontainer" style="font-size:1.2rem;">

      <div id="choose-data" class="pt-4">
				<div id="choose-tissue"></div>
				<div id="choose-blood"></div>
        <hr>
      </div>
    
      <div id="display-meta">

        <div id="" class="row" style="min-height: min(480px, 40vw);">

            <div class="col-md-6 col-lg-6">
              <div class="row">
                <div class="col-12"><h4>Metadata</h4></div>
                <div id="embedding-meta" class="col-12"></div>
              </div>
              <div id="meta-controls" class="row mt-2"></div>
            </div>

            <div class="col-md-6 col-lg-6">
              <div class="row">
                <div class="col-12"><h4>Gene Expression</h4></div>
                <div id="gene-loading-spinner" class="col-2 offset-5"><div class="loader"></div></div>
                <div id="embedding-gene" class="col-12"></div>
              </div>
              <div id="gene-controls" class="row mt-2"></div>
            </div>

            <div class="col-md-4 col-lg-4">
              <div class="row">
        <div class="col-12"><h4>Expression Summary</h4></div>
                <div id="gene-bars" class="col"></div>
                <div class="w-100"></div>
                <div id="gene-legend" class="col"></div>
              </div>
            </div>

            <div class="col-md-8 col-lg-8">
              <div class="row">
                <div class="col">
                  <h4>Gene Table</h4>
                  <div id="mytable"></div>
                </div>
              </div>
            </div>

        </div> <!-- row -->

      </div> <!-- display-meta -->

          <!--
            <div>
              <h3>Heatmap of cell counts</h3>
              <a class="btn btn-secondary" data-toggle="collapse" href="#collapse-heatmap" role="button" aria-expanded="false" aria-controls="collapse-heatmap">Show / Hide</a>
            </div>
            <div class="collapse show" id="collapse-heatmap">
              <div id="meta-row-2" class="row">
                <div id="meta-heatmap" class="col-12"></div>
              </div>
            </div>
          -->

      <div class="row">
          <div id="geneinfo" class="col-12"></div>
      </div>

      <div id="cluster-abundance">
        <div class="row">
          <div class="col"><h4>Abundance of Cell Clusters</h4></div>
        </div>
        <div class="row pb-5">
          <div class="col-2"></div>
          <div id="meta-boxplot" class="col-6"></div>
          <div class="col-2"></div>
        </div>
      </div>

      <div class="row my-4">
        <div class="col">
          <h2 class="f2 my-3">Cite this data</h2>
          <p class="lh-copy">Thomas MF, Slowikowski K, Manakongtreecheep K, Sen P, Tantivit J, Nasrallah M, et al. <a class="link blue" href="https://doi.org/10.1101/2021.09.17.460868">Altered interactions between circulating and tissue-resident CD8 T cells with the colonic mucosa define colitis associated with immune checkpoint inhibitors.</a> bioRxiv. 2021. p. 2021.09.17.460868. doi:10.1101/2021.09.17.460868</p>
          <p>
          <details><summary>BibTex</summary><pre><code>
  @UNPUBLISHED{Thomas2021,
    title    = "{Altered interactions between circulating and tissue-resident CD8
                T cells with the colonic mucosa define colitis associated with
                immune checkpoint inhibitors}",
    author   = "Thomas, Molly Fisher and Slowikowski, Kamil and
                Manakongtreecheep, Kasidet and Sen, Pritha and Tantivit, Jessica
                and Nasrallah, Mazen and Smith, Neal P and Ramesh, Swetha and
                Zubiri, Leyre and Tirard, Alice and Arnold, Benjamin Y and
                Nieman, Linda T and Chen, Jonathan H and Eisenhaure, Thomas and
                Pelka, Karin and Xu, Katherine H and Jorgji, Vjola and Pinto,
                Christopher J and Sharova, Tatyana and Glasser, Rachel and Chan,
                Elaina Puiyee and Sullivan, Ryan J and Khalili, Hamed and Juric,
                Dejan and Boland, Genevieve M and Dougan, Michael and Hacohen,
                Nir and Reynolds, Kerry L and Li, Bo and Villani,
                Alexandra-Chlo{\'e}",
    journal  = "bioRxiv",
    pages    = "2021.09.17.460868",
    month    =  sep,
    year     =  2021,
    language = "en",
    doi      = "10.1101/2021.09.17.460868"
  }
            </code></pre></details>
          </p>
        </div>
      </div>

    </div>

    <div class="bg-light py-5">
      <footer class="text-center mastfoot my-5">
        <div id="contact-me" class="inner">
          <p><a href="https://cell.guide">Cell Guide</a> is a project by <a href="https://twitter.com/slowkow">@slowkow</a>.</p>
        </div>
      </footer>
    </div>
    `
  )

  var renderer = function() {

    function drawMetaHex(fieldName, onDone = () => null) {
      console.log(`drawMetaHex(fieldName = ${fieldName})`);

      //var fieldName = 'cluster';
      var metaInfo = db.findMetaInfo(fieldName);

      if (metaInfo && metaInfo.type == "enum") {
        var meta_valcounts = get_groups(fieldName)
        var this_pal = metaInfo.valCounts.length <= pals.okabe.length ?
          pals.okabe : pals.mpn65
        if (fieldName === 'case' && metaInfo.valCounts.length === 2) {
          this_pal = pals.okabe.slice(0,2).reverse()
        }
        var metaColors = {}
        // {'healthy': {r: 240, g: 160, b: 255, opacity: 1}, ...}
        for (var i = 0; i < meta_valcounts.length; i++) {
          var x = meta_valcounts[i][0]
          var color_i = i
          if (color_i > this_pal.length - 1) {
            color_i = color_i % this_pal.length
          }
          metaColors[x] = d3.color(this_pal[color_i])
        }
      }

      // _dump(metaColors)

      function doDrawMetaHex() {
        //
        let data = g_mydata
        //
        let canvas_scale = 2
        let radius = g_size.hex_radius
        let legend_width = 90
        let plot_width = 360 * canvas_scale
        let plot_height = 360 * canvas_scale
        //
        let legend_margin = {
          top: 35, right: 0, bottom: 0, left: 5
        }
        let margin = {
          top: 35, right: 10, bottom: 10, left: 10
        }
        //
        const container = d3.select("#embedding-meta")
        container.html("")
        const canvas_id = "mycanvas-meta-plot"
        // d3.select(canvas_id).remove()
        const canvas = container.append("canvas")
          .attr("id", canvas_id)
          .style("width", "100%")
          .node()
        canvas.width = plot_width + legend_width * canvas_scale
        canvas.height = plot_height

        const context = canvas.getContext('2d')

        const svg = container.append('svg')
          .attr("viewBox",
            `0 0 ${canvas.width / canvas_scale} ${1 + canvas.height / canvas_scale}`
          )
          .style("position", "absolute")
          .style("left", "10px")
          .style("top", "0px")
          .style("height", "100%")

        let meta_title = metaInfo.valCounts ?
          `${metaInfo.label} (n = ${metaInfo.valCounts.length})` :
          `${metaInfo.label}`

        // const font_size = 10
        let font_size = g_font_size

        svg.append("text")
          .attr("x", margin.left)
          .attr("y", margin.top - 18)
          .attr("font-family", "sans-serif")
          .attr("font-size", `${font_size * 1.5}px`)
          .text(meta_title)
        
        svg.append("text")
          .attr("x", margin.left / 2)
          .attr("y", plot_height / canvas_scale - margin.bottom / 2)
          .attr("font-family", "sans-serif")
          .attr("font-size", `${font_size * 1.4}px`)
          .text(d3.format(",")(db.conf.sampleCount) + " cells")

        let panelBorder = svg.append("rect")
          .attr("x", margin.left - 9)
          .attr("y", margin.top - 10)
          .attr("width", plot_width / canvas_scale - margin.right / 2)
          .attr("height", plot_height / canvas_scale - margin.top + margin.bottom)
          .style("stroke", "black")
          .style("fill", "none")
          .style("stroke-width", "0.5px");

      var x = d3.scaleLinear()
        .domain(d3.extent(data, d => d.x))
        .range([margin.left * canvas_scale, plot_width - margin.right * canvas_scale])

      var y = d3.scaleLinear()
        .domain(d3.extent(data, d => d.y))
        .range([-10 + plot_height - margin.bottom * canvas_scale, margin.top * canvas_scale])

      var hexbin = d3.hexbin()
        .x(d => x(d.x))
        .y(d => y(d.y))
        .radius(radius * (plot_width - legend_width * canvas_scale) / (plot_height))
        .extent([
          [margin.left, margin.top],
          [plot_width - margin.right, plot_height - margin.bottom]
        ])
        var bins = hexbin(g_mydata)

        // var color = d3.scaleSequential(d3.interpolateBuPu)
        //   .domain([0, d3.max(bins, bin => bin.length) / 2])
        
        context.fillStyle = "#fff";

        // context.strokeStyle = "black";
        // context.strokeRect(
        //   margin.left - 10,
        //   margin.top - 10,
        //   plot_width,
        //   plot_height - margin.top + margin.bottom
        // )

        var hex = new Path2D(hexbin.hexagon())

        if (metaInfo.type == "enum") {

          // mixing colors...
          // TODO: consider using HCL space https://bl.ocks.org/mbostock/3014589
          function colorEnum(bin) {
            // Input: [{health: "healthy"}, {health: "not healthy"}, ...]
            // Ouput: [["healthy", 0.25], ["not healthy", 0.75]]
            var a = Array.from(d3.rollup(
              bin, v => v.length / bin.length, d => d[fieldName]
            ))
            var r = 0, g = 0, b = 0
            for (var i = 0; i < a.length; i++) {
              r += a[i][1] * metaColors[a[i][0]].r
              g += a[i][1] * metaColors[a[i][0]].g
              b += a[i][1] * metaColors[a[i][0]].b
            }
            return(d3.rgb(r, g, b).toString())
          }

          bins.forEach(function(bin) {
            context.translate(bin.x, bin.y)
            context.fillStyle = colorEnum(bin)
            context.fill(hex)
            context.strokeStyle = context.fillStyle
            context.lineWidth = 0
            context.stroke(hex)
            context.setTransform(1, 0, 0, 1, 0, 0)
          });

          var ordinal = d3.scaleOrdinal()
            // .domain(meta_valcounts.map(d => d[0]))
            .domain(meta_valcounts.map(d => `${d[0]} (n = ${d3.format(",d")(d[1])})`))
            .range(meta_valcounts.map(d => metaColors[d[0]]))

          // meta_valcounts.map(d => [`${d[0]} (n = ${d3.format(",d")(d[1])})`, d[1]])

          svg.append("g")
            .attr("class", "legendOrdinal")
            .attr("transform",
              `translate(
                ${plot_width / canvas_scale - margin.right / canvas_scale + legend_margin.left},
                ${legend_margin.top})`
            )

          var legendOrdinal = svg => {
            const box_height = Math.min(
              16,
              (0.9 * plot_height / canvas_scale) / meta_valcounts.length
            )
            let font_size = Math.min(g_font_size, box_height)
            const g = svg
                .attr("text-anchor", "start")
                .attr("font-family", "sans-serif")
                // .attr("font-size", 14)
                .attr("font-size", `${font_size}px`)
              .selectAll("g")
              .data(ordinal.domain())
              .join("g")
                .attr("transform", (d, i) => `translate(0,${i * box_height})`)
            g.append("rect")
              .attr("width", 10)
              .attr("height", box_height * 0.8)
              .attr("fill", ordinal)
            g.append("text")
              .attr("x", 10 + 3)
              .attr("y", box_height / 2)
              .attr("alignment-baseline", "middle")
              // .attr("dy", "0.35em")
              .text(d => d)
          }

          // var legendOrdinal = d3.legendColor()
          //   .shape("path", d3.symbol().type(d3.symbolCircle).size(200)())
          //   // .labelWrap(legend_width - 30)
          //   // .shapePadding(10)
          //   .scale(ordinal);

          // if (metaColors.length < 10) {
            svg.select(".legendOrdinal")
              .call(legendOrdinal);
          // }
          
          // function text_position(d,i) {
          //   var c = 16;   // number of rows
          //   var h = 20;  // legend entry height
          //   var w = 150; // legend entry width (so we can position the next column) 
          //   var tx = 10; // tx/ty are essentially margin values
          //   var ty = 10;
          //   var y = i % c * h + ty;
          //   var x = Math.floor(i / c) * w + tx;
          //   return "translate(" + x + "," + y + ")";
          // }

          // svg.select(".legendOrdinal")
          //   .attr("font-family", "sans-serif")
          //   .attr("font-size", `${plot_height / metaColors.length}px`)
          //   .append("g")
          //   .selectAll("g")
          //   .data(items)
          //   .join("g")
          //   .append("text")
          //   // .attr("transform", (d, i) => `translate(0,${i * 16})`)
          //   .attr("transform", text_position)
          //   .text(d => d[0])


        } else if (metaInfo.type == "int" || metaInfo.type == "float") {

          const color_domain = [
            0, d3.max(bins, bin => d3.mean(bin, d => d[fieldName]))
          ]

          // var color = d3.scaleSequential(d3.interpolateBuPu)
          //   .domain(color_domain)

          // var color = d3.scaleQuantize()
          //   .domain(color_domain)
          //   .range([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0].map(d => d3.interpolateBuPu(d)))

          // let color_range = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
          
          // // Quantile scale (weird?)
          // let color_range = linspace(0, 1, 11)
          // var color = d3.scaleQuantile()
          //   .domain(bins.map(bin => d3.mean(bin, d => d[fieldName])))
          //   .range(color_range.map(d => d3.interpolateBuPu(d)))


          // Sequential scale (works)
          var extent = d3.extent(bins, bin => d3.mean(bin, d => d[fieldName]))
          var color = null
          var format = null
          if (extent[0] < 0 && extent[1] > 0) {
            color = d3.scaleDiverging(t => d3.interpolateRdBu(1 - t))
              .domain([extent[0], 0, extent[1]])

          } else {
            color = d3.scaleSequential(d3.interpolateBuPu)
              .domain([extent[0], extent[1]])
          }
          var color_range = linspace(extent[0], extent[1], 10)

          bins.forEach(function(bin) {
            context.translate(bin.x, bin.y)
            context.fillStyle = color(d3.mean(bin, d => d[fieldName]))
            context.fill(hex)
            context.setTransform(1, 0, 0, 1, 0, 0)
          });

          // svg.append("text")
          //   .attr("x", plot_width / canvas_scale + legend_margin.left)
          //   .attr("y", legend_margin.top)
          //   .attr("font-family", "sans-serif")
          //   .attr("font-size", `${g_font_size}px`)
          //   .attr("font-weight", "bold")
          //   .text(metaInfo.label)

          vertical_legend({
            svg: svg,
            color: color,
            title: metaInfo.label,
            offsetLeft: plot_width / canvas_scale - margin.right / canvas_scale + legend_margin.left,
            offsetTop: legend_margin.top
          })

          //var legend = svg.selectAll("g.legend_colorbar")
          //  // .data(color.range().reverse())
          //  .data(color_range)
          //  .enter()
          //  .append("g")
          //  .attr("class", "legend_colorbar")
          ////
          //let legend_rect_height = (
          //  plot_height / canvas_scale - legend_margin.top - legend_margin.bottom
          //) / (color_range.length + 2)
          ////
          //legend
          //  .append('rect')
          //  .attr("x", plot_width / canvas_scale + legend_margin.left)
          //  .attr("y", function(d, i) {
          //    return legend_margin.top + (color_range.length - i - 0.5) * legend_rect_height
          //  })
          // .attr("width", 15)
          // .attr("height", legend_rect_height)
          // .style("fill", d => color(d))
          ////
          //legend
          //  .append('text')
          //  .attr("x", plot_width / canvas_scale + legend_margin.left + 20)
          //  .attr("y", function(d, i) {
          //    return legend_margin.top + (color_range.length - i) * legend_rect_height
          //  })
          //  .attr("alignment-baseline", "middle")
          //  .style("font-size", `${g_font_size}px`)
          //  .text(d => d3.format(",.2r")(d))

        }

        // cluster labels

        var clusterField = db.conf.clusterField
        let clusterLabels = d3.nest()
          .key(d => d[clusterField])
          .rollup(v => [d3.mean(v, d => d.x), d3.mean(v, d => d.y)])
          .entries(g_mydata)

        // var x2 = d3.scaleLinear()
        //   .domain(d3.extent(g_mydata, d => d.x))
        //   .range([margin.left * 10, plot_width])
        var x2 = d3.scaleLinear()
          .domain(d3.extent(g_mydata, d => d.x))
          .range([margin.left, plot_width - legend_width / canvas_scale])

        var y2 = d3.scaleLinear()
          .domain(d3.extent(g_mydata, d => d.y))
          .rangeRound([plot_height - margin.bottom, margin.top * canvas_scale])

        function x_cl(d) {
          return (x2(d.value[0]) - margin.right) / canvas_scale
        }

        svg.selectAll("g.cluster_label")
          .data(clusterLabels)
          .enter()
          .append("g")
          .append("text")
          .attr("class","cluster_label")
          .attr("x", d => x2(d.value[0]) / canvas_scale)
          .attr("y", d => y2(d.value[1]) / canvas_scale)
          .style("display", $("#show-cluster-labels").is(":checked") ? "block" : "none")
          .text(d => d.key)

        drawBars(metaInfo)

      } // doDrawMetaHex

      function doDrawMetaHex1() {
        // TODO The xy coordinates in this function need a total rewrite.
        var radius = g_size.hex_radius
        var legend_width = 280
        const canvas_scale = 2
        var plot_width = 520 * canvas_scale
        var plot_height = 460 * canvas_scale
        var legend_margin = {
          top: 40, right: 0, bottom: 0, left: 5
        }
        var margin = {
          top: 45, right: 90, bottom: 10, left: 0
        }

        d3.select("#embedding-meta").html("")
          // .style("position", "relative")
          // .style("display", "inline-block")
        // const container = d3.select("#mydiv-meta").append("div")
        //   .attr("id", "embedding-meta")
        //   .style("position", "relative")
        //   .style("display", "inline-block")
        //   // .style("display", "none")
        // const canvas = container.append('canvas').node()
        const container = d3.select("#embedding-meta")
        const canvas = container.append("canvas")
          .attr("id", "mycanvas-meta-plot")
          .style("width", "100%")
          .style("margin-left", '10px')
          // .style("width", `${plot_width + legend_width}px`)
          .node()
        const context = canvas.getContext('2d')

        canvas.width = plot_width + legend_width * canvas_scale
        canvas.height = plot_height * (canvas_scale * 0.6)

        const svg_width = canvas.width / canvas_scale - legend_width / 1.7 
        const svg_height = canvas.height / canvas_scale
        const svg = container.append('svg')
          // .attr("width", svg_width)
          // .attr("height", svg_height)
          .attr("viewBox", `-10 0 ${svg_width} ${svg_height}`)
          .style("position", "absolute")
          .style("top", '0px')
          .style("left", '0px')
          // .style("left", `${margin.left}px`)

        let meta_title = metaInfo.valCounts ?
          `${metaInfo.label} (n = ${metaInfo.valCounts.length})` :
          `${metaInfo.label}`

        svg.append("text")
          .attr("x", margin.left)
          .attr("y", margin.top - 20)
          .attr("font-family", "sans-serif")
          .attr("font-size", `${g_font_size}px`)
          .text(meta_title)

        svg.append("text")
          .attr("x", margin.left / 2)
          .attr("y", plot_height / canvas_scale - margin.bottom)
          .attr("font-family", "sans-serif")
          .attr("font-size", `${g_font_size}px`)
          .text(d3.format(",")(db.conf.sampleCount) + " cells")

        let panelBorder = svg.append("rect")
          .attr("x", margin.left - 10)
          .attr("y", margin.top - 10)
          .attr("width", plot_width / canvas_scale - margin.right / canvas_scale)
          .attr("height", plot_height / canvas_scale - margin.top + margin.bottom)
          .style("stroke", "black")
          .style("fill", "none")
          .style("stroke-width", "0.5px");

        var x = d3.scaleLinear()
          .domain(d3.extent(g_mydata, d => d.x))
          .range([margin.left, plot_width + margin.right])

        var y = d3.scaleLinear()
          .domain(d3.extent(g_mydata, d => d.y))
          .rangeRound([plot_height - margin.bottom, margin.top])

        var hexbin = d3.hexbin()
          .x(d => x(d.x))
          .y(d => y(d.y))
          .radius(radius * (plot_width - legend_width) / plot_height + 0.2)
          .extent([[margin.left, margin.top], [plot_width - margin.right, plot_height - margin.bottom]])

        var bins = hexbin(g_mydata)

        // var color = d3.scaleSequential(d3.interpolateBuPu)
        //   .domain([0, d3.max(bins, bin => bin.length) / 2])
        
        context.fillStyle = "#fff";

        // context.strokeStyle = "black";
        // context.strokeRect(
        //   margin.left - 10,
        //   margin.top - 10,
        //   plot_width,
        //   plot_height - margin.top + margin.bottom
        // )

        var hex = new Path2D(hexbin.hexagon());

        if (metaInfo.type == "enum") {

          // mixing colors...
          // TODO: consider using HCL space https://bl.ocks.org/mbostock/3014589
          function colorEnum(bin) {
            // Input: [{health: "healthy"}, {health: "not healthy"}, ...]
            // Ouput: [["healthy", 0.25], ["not healthy", 0.75]]
            var a = Array.from(d3.rollup(
              bin, v => v.length / bin.length, d => d[fieldName]
            ))
            var r = 0, g = 0, b = 0
            for (var i = 0; i < a.length; i++) {
              r += a[i][1] * metaColors[a[i][0]].r
              g += a[i][1] * metaColors[a[i][0]].g
              b += a[i][1] * metaColors[a[i][0]].b
            }
            return(d3.rgb(r, g, b).toString())
          }

          bins.forEach(function(bin) {
            context.translate(bin.x, bin.y)
            context.fillStyle = colorEnum(bin)
            context.fill(hex)
            context.strokeStyle = context.fillStyle
            context.lineWidth = 0
            context.stroke(hex)
            context.setTransform(1, 0, 0, 1, 0, 0)
          });

          var ordinal = d3.scaleOrdinal()
            // .domain(meta_valcounts.map(d => d[0]))
            .domain(meta_valcounts.map(d => `${d[0]} (n = ${d3.format(",d")(d[1])})`))
            .range(meta_valcounts.map(d => metaColors[d[0]]))

          // meta_valcounts.map(d => [`${d[0]} (n = ${d3.format(",d")(d[1])})`, d[1]])

          svg.append("g")
            .attr("class", "legendOrdinal")
            .attr("transform",
              `translate(
                ${plot_width / canvas_scale - margin.right / canvas_scale + legend_margin.left},
                ${legend_margin.top})`
            )

          var legendOrdinal = svg => {
            const box_height = Math.min(
              20,
              (0.9 * plot_height / canvas_scale) / meta_valcounts.length
            )
            let font_size = Math.min(g_font_size, box_height)
            const g = svg
                .attr("text-anchor", "start")
                .attr("font-family", "sans-serif")
                // .attr("font-size", 14)
                .attr("font-size", `${font_size}px`)
              .selectAll("g")
              .data(ordinal.domain())
              .join("g")
                .attr("transform", (d, i) => `translate(0,${i * box_height})`)
            g.append("rect")
              .attr("width", 10)
              .attr("height", box_height * 0.9)
              .attr("fill", ordinal)
            g.append("text")
              .attr("x", 10 + 3)
              .attr("y", box_height / 2)
              .attr("alignment-baseline", "middle")
              // .attr("dy", "0.35em")
              .text(d => d)
          }

          // var legendOrdinal = d3.legendColor()
          //   .shape("path", d3.symbol().type(d3.symbolCircle).size(200)())
          //   // .labelWrap(legend_width - 30)
          //   // .shapePadding(10)
          //   .scale(ordinal);

          // if (metaColors.length < 10) {
            svg.select(".legendOrdinal")
              .call(legendOrdinal);
          // }
          
          // function text_position(d,i) {
          //   var c = 16;   // number of rows
          //   var h = 20;  // legend entry height
          //   var w = 150; // legend entry width (so we can position the next column) 
          //   var tx = 10; // tx/ty are essentially margin values
          //   var ty = 10;
          //   var y = i % c * h + ty;
          //   var x = Math.floor(i / c) * w + tx;
          //   return "translate(" + x + "," + y + ")";
          // }

          // svg.select(".legendOrdinal")
          //   .attr("font-family", "sans-serif")
          //   .attr("font-size", `${plot_height / metaColors.length}px`)
          //   .append("g")
          //   .selectAll("g")
          //   .data(items)
          //   .join("g")
          //   .append("text")
          //   // .attr("transform", (d, i) => `translate(0,${i * 16})`)
          //   .attr("transform", text_position)
          //   .text(d => d[0])


        } else if (metaInfo.type == "int" || metaInfo.type == "float") {

          const color_domain = [
            0, d3.max(bins, bin => d3.mean(bin, d => d[fieldName]))
          ]

          // var color = d3.scaleSequential(d3.interpolateBuPu)
          //   .domain(color_domain)

          // var color = d3.scaleQuantize()
          //   .domain(color_domain)
          //   .range([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0].map(d => d3.interpolateBuPu(d)))

          // let color_range = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
          
          // // Quantile scale (weird?)
          // let color_range = linspace(0, 1, 11)
          // var color = d3.scaleQuantile()
          //   .domain(bins.map(bin => d3.mean(bin, d => d[fieldName])))
          //   .range(color_range.map(d => d3.interpolateBuPu(d)))


          // Sequential scale (works)
          var extent = d3.extent(bins, bin => d3.mean(bin, d => d[fieldName]))
          var color = null
          var format = null
          if (extent[0] < 0 && extent[1] > 0) {
            color = d3.scaleDiverging()
              .domain([extent[0], 0, extent[1]])
              // .interpolator(d3.interpolatePiYG)
              .interpolator(d3.interpolateRdYlBu)
          } else {
            color = d3.scaleSequential(d3.interpolateBuPu)
              .domain([extent[0], extent[1]])
          }
          var color_range = linspace(extent[0], extent[1], 10)

          bins.forEach(function(bin) {
            context.translate(bin.x, bin.y)
            context.fillStyle = color(d3.mean(bin, d => d[fieldName]))
            context.fill(hex)
            context.setTransform(1, 0, 0, 1, 0, 0)
          });

          // svg.append("text")
          //   .attr("x", plot_width / canvas_scale + legend_margin.left)
          //   .attr("y", legend_margin.top)
          //   .attr("font-family", "sans-serif")
          //   .attr("font-size", `${g_font_size}px`)
          //   .attr("font-weight", "bold")
          //   .text(metaInfo.label)

          vertical_legend({
            svg: svg,
            color: color,
            title: metaInfo.label,
            offsetLeft: plot_width / canvas_scale - margin.right / canvas_scale + legend_margin.left,
            offsetTop: legend_margin.top
          })

          //var legend = svg.selectAll("g.legend_colorbar")
          //  // .data(color.range().reverse())
          //  .data(color_range)
          //  .enter()
          //  .append("g")
          //  .attr("class", "legend_colorbar")
          ////
          //let legend_rect_height = (
          //  plot_height / canvas_scale - legend_margin.top - legend_margin.bottom
          //) / (color_range.length + 2)
          ////
          //legend
          //  .append('rect')
          //  .attr("x", plot_width / canvas_scale + legend_margin.left)
          //  .attr("y", function(d, i) {
          //    return legend_margin.top + (color_range.length - i - 0.5) * legend_rect_height
          //  })
          // .attr("width", 15)
          // .attr("height", legend_rect_height)
          // .style("fill", d => color(d))
          ////
          //legend
          //  .append('text')
          //  .attr("x", plot_width / canvas_scale + legend_margin.left + 20)
          //  .attr("y", function(d, i) {
          //    return legend_margin.top + (color_range.length - i) * legend_rect_height
          //  })
          //  .attr("alignment-baseline", "middle")
          //  .style("font-size", `${g_font_size}px`)
          //  .text(d => d3.format(",.2r")(d))

        }

        // cluster labels

        var clusterField = db.conf.clusterField
        let clusterLabels = d3.nest()
          .key(d => d[clusterField])
          .rollup(v => [d3.mean(v, d => d.x), d3.mean(v, d => d.y)])
          .entries(g_mydata)

        let cluster_label_size = "1em"

        var x2 = d3.scaleLinear()
          .domain(d3.extent(g_mydata, d => d.x))
          .range([margin.left * 10, plot_width])

        var y2 = d3.scaleLinear()
          .domain(d3.extent(g_mydata, d => d.y))
          .rangeRound([plot_height - margin.bottom, margin.top * canvas_scale])

        function x_cl(d) {
          return (x2(d.value[0]) - margin.right) / canvas_scale
        }

        svg.selectAll("g.cluster_label")
          .data(clusterLabels)
          .enter()
          .append("g")
          .append("text")
          .attr("class","cluster_label")
          .attr("x", x_cl)
          .attr("y", d => y2(d.value[1]) / canvas_scale)
          // .attr("text-anchor", "middle")
          // .attr("font-family", "sans-serif")
          // .attr("font-size", cluster_label_size)
          // // .attr("font-weight", "bold")
          // .style("fill", "black")
          // .style("stroke", "white")
          // .style("stroke-width", "4px")
          // .style("stroke-align", "outset")
          .style("display", $("#show-cluster-labels").is(":checked") ? "block" : "none")
          .text(d => d.key)

        // svg.selectAll("g.cluster_label")
        //   .data(clusterLabels)
        //   .enter()
        //   .append("g")
        //   .append("text")
        //   .attr("class","cluster_label")
        //   .attr("x", x_cl)
        //   .attr("y", d => y2(d.value[1]) / canvas_scale)
        //   .attr("text-anchor", "middle")
        //   .attr("font-family", "sans-serif")
        //   .attr("font-size", cluster_label_size)
        //   .attr("font-weight", "bold")
        //   .style("fill", "black")
        //   .style("display", $("#show-cluster-labels").is(":checked") ? "block" : "none")
        //   .text(d => d.key)

        drawBars(metaInfo)

      } // doDrawMetaHex

      function gotMetaArr(metaArr, metaInfo) {
        var fieldName = metaInfo.name;
        if (metaInfo.valCounts) {
          for (var i = 0; i < metaArr.length; i++) {
            g_mydata[i][fieldName] = metaInfo.valCounts[metaArr[i]][0]
          }
        } else if (metaInfo.origVals){
          for (var i = 0; i < metaInfo.origVals.length; i++) {
            g_mydata[i][fieldName] = metaInfo.origVals[i]
          }
        }
        doDrawMetaHex()
        onDone()
      }

      if (metaInfo) {
        if (metaInfo.origVals) {
          // for numeric fields, the raw data is already in memory
          gotMetaArr(metaInfo.origVals, metaInfo)
        } else {
          // other fields may not be loaded yet
          db.loadMetaVec(metaInfo, gotMetaArr, onProgressConsole);
        }
      }

      // db.loadMetaVec(metaInfo, gotMetaArr, onProgressConsole);
    } // drawMetaHex()

    function draw_gene_bars(geneSym, gene_groupby) {

      var metaInfo = null

      if (gene_groupby !== "none") {
        metaInfo = db.findMetaInfo(gene_groupby)
        if (metaInfo.type != "enum" || metaInfo.valCounts.length > pals.okabe.length) {
          return
        }
      } else {
        return draw_gene_bars_none(geneSym)
      }

      let groupKey = db.conf.clusterField
      let groups = get_groups(groupKey)
      let group_names = groups.map(d => d[0])

      let subgroupKey = metaInfo.name
      let subgroupLevels = metaInfo.valCounts.map(d => d[0]).slice().sort()

      var data = []
      var x_min = Infinity
      var x_max = -Infinity
      for (var i = 0; i < groups.length; i++) {
        var res = {}
        res[groupKey] = groups[i][0]
        var counts = Array.from(d3.rollup(
          g_mydata.filter(d => d[groupKey] == groups[i][0]),
          v => v.filter(d => d.gene > 0).length / v.length, d => d[subgroupKey]
        ))
        for (var j = 0; j < counts.length; j++) {
          let val = counts[j][1]
          res[counts[j][0]] = val
          x_min = val < x_min ? val : x_min
          x_max = val > x_max ? val : x_max
        }
        data.push(res)
      }

      var el = $("#gene-bars")
      el.html("")
      const width = g_size.gene_bars.width
      const height = g_size.gene_bars.height(groups.length, subgroupLevels.length)
      const longest_group_key = d3.greatest(
        group_names, (a, b) => d3.ascending(a.length, b.length)
      )
      const margin = {
        top: 15, right: 15, bottom: 40,
        left: getTextWidth(longest_group_key, "15px arial") + 9
      }

      let y0 = d3.scaleBand()
        .domain(data.map(d => d[groupKey]))
        // .rangeRound([margin.top, height - margin.bottom])
        .range([margin.top + 2, height - margin.bottom - 2])
        .paddingInner(0.1)
        .paddingOuter(0.0)

      let y1 = d3.scaleBand()
        .domain(subgroupLevels)
        .rangeRound([0, y0.bandwidth()])
        .padding(0.05)

      let x = d3.scaleLinear()
        // .domain([0, d3.max(data, d => d3.max(subgroupLevels, key => d[key]))]).nice()
        .domain([0, 1]).nice()
        .rangeRound([margin.left + 5, width - margin.right - 5])

      let color_range = subgroupLevels.length <= 8 ? pals.okabe : pals.mpn65
      if (gene_groupby === 'case') {
        color_range = pals.okabe.slice(0,2).reverse()
      }
      let color = d3.scaleOrdinal()
        .domain(subgroupLevels)
        .range(color_range)

      let xAxis = g => g
        .attr("transform", `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).tickSizeOuter(0).ticks(2, "~p"))
        .call(g => g.select(".domain").remove())
        .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))
        // .call(
        //   g => g.select(".tick:last-of-type text").clone()
        //     .attr("x", 15)
        //     .attr("text-anchor", "start")
        //     // .attr("font-weight", "bold")
        //     .attr("font-size", `${g_font_size}px`)
        //     .text("cells")
        // )

      let yAxis = g => g
        .attr("transform", `translate(${margin.left},0)`)
        .call(d3.axisLeft(y0).ticks(null, "s"))
        .call(g => g.select(".domain").remove())
        .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))

      const legend = svg => {
        const g = svg
            .attr("transform", `translate(${width},0)`)
            .attr("text-anchor", "end")
            .attr("font-family", "sans-serif")
            .attr("font-size", 14)
          .selectAll("g")
          .data(color.domain())
          .join("g")
            .attr("transform", (d, i) => `translate(0,${i * 20})`);
        g.append("rect")
            .attr("x", -19)
            .attr("width", 19)
            .attr("height", 19)
            .attr("fill", color);
        g.append("text")
            .attr("x", -24)
            .attr("y", 9.5)
            .attr("dy", "0.35em")
            .text(d => d);
      }

      const svg = d3.select("#gene-bars").append('svg')
        .attr("id", "gene-bars-svg")
        .attr("viewBox", `0 0 ${g_size.gene_bars.width + g_size.gene_boxplot.width} ${height}`)

      // background stripes
      const bg = svg.append("g")
        .selectAll("g")
        .data(data)
        .join("g")
      bg.append("rect")
        .attr("fill", (d, i) => i % 2 == 0 ? "#ffffff00" : "#eee")
        .attr("x", d => x(0) - 4)
        .attr("y", d => y0(d[groupKey]) - y1.bandwidth() / 8)
        .attr("height", y0.bandwidth() + y1.bandwidth() / 4)
        .attr("width", d => x(1) - x(0) + 8)

      let panelBorder = svg.append("rect")
        .attr("x", margin.left)
        .attr("y", margin.top)
        .attr("height", height - margin.top - margin.bottom)
        .attr("width", width - margin.left - margin.right)
        .style("stroke", "black")
        .style("fill", "none")
        .style("stroke-width", "0.5px");

      svg.append("g")
        .selectAll("g")
        .data(data)
        .join("g")
          .attr("transform", d => `translate(0,${y0(d[groupKey])})`)
        .selectAll("rect")
        .data(d => subgroupLevels.map(key => ({key, value: d[key]})))
        .join("rect")
          .attr("x", d => x(0))
          .attr("y", d => y1(d.key))
          .attr("height", y1.bandwidth())
          .attr("width", d => x(d.value) - x(0))
          .attr("fill", d => color(d.key))
          .attr("stroke", "#000000")
          .attr("stroke-width", "0.5px")

      svg.append("text")
        .attr("x", margin.left)
        .attr("y", margin.top - 5)
        .attr("font-family", "sans-serif")
        .attr("font-size", `${g_font_size}px`)
        // .attr("font-style", "italic")
        // .html(`Cells with <tspan font-style="italic">${geneSym}</tspan>`)
        .html(`Cells`)

      svg.append("g")
          .call(xAxis);

      svg.append("g")
          .call(yAxis);

      // svg.append("g")
      //     .call(legend);
    }

    function draw_gene_bars_none(geneSym) {

      let groupKey = db.conf.clusterField
      let groups = get_groups(groupKey)
      let group_names = groups.map(d => d[0])

      // console.log(`draw_gene_bars_none('${geneSym}') groupKey = ${groupKey}`)

      var data = []
      var x_min = Infinity
      var x_max = -Infinity
      for (var i = 0; i < groups.length; i++) {
        let items = g_mydata.filter(d => d[groupKey] == groups[i][0])
        let val = items.filter(d => d.gene > 0).length / items.length
        x_min = val < x_min ? val : x_min
        x_max = val > x_max ? val : x_max
        let datum = {}
        datum[groupKey] = groups[i][0]
        datum.value = val
        data.push(datum)
      }

      var el = $("#gene-bars")
      el.html("")
      const width = g_size.gene_bars.width
      const height = g_size.gene_bars.height(groups.length, 1)
      const longest_group_key = d3.greatest(
        group_names, (a, b) => d3.ascending(a.length, b.length)
      )
      let margin = {
        top: 15, right: 15, bottom: 40,
        left: getTextWidth(longest_group_key, "15px arial") + 9
      }

      let y0 = d3.scaleBand()
        .domain(data.map(d => d[groupKey]))
        // .rangeRound([margin.top, height - margin.bottom])
        .range([margin.top + 2, height - margin.bottom - 2])
        .paddingInner(0.1)
        .paddingOuter(0.0)

      let x = d3.scaleLinear()
        .domain([0, 1]).nice()
        .rangeRound([margin.left + 5, width - margin.right - 5])

      let xAxis = g => g
        .attr("transform", `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).tickSizeOuter(0).ticks(2, "~p"))
        .call(g => g.select(".domain").remove())
        .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))

      let yAxis = g => g
        .attr("transform", `translate(${margin.left},0)`)
        .call(d3.axisLeft(y0).ticks(null, "s"))
        .call(g => g.select(".domain").remove())
        .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))

      const svg = d3.select("#gene-bars").append('svg')
        .attr("viewBox", `0 0 ${g_size.gene_bars.width + g_size.gene_boxplot.width} ${height}`)
        .attr("id", "gene-bars-svg")

      // background stripes
      const bg = svg.append("g")
        .selectAll("g")
        .data(data)
        .join("g")
      bg.append("rect")
        .attr("fill", (d, i) => i % 2 == 0 ? "#ffffff00" : "#eee")
        .attr("x", d => x(0) - 4)
        .attr("y", d => y0(d[groupKey]))
        .attr("height", y0.bandwidth())
        .attr("width", d => x(1) - x(0) + 8)

      let panelBorder = svg.append("rect")
        .attr("x", margin.left)
        .attr("y", margin.top)
        .attr("height", height - margin.top - margin.bottom)
        .attr("width", width - margin.left - margin.right)
        .style("stroke", "black")
        .style("fill", "none")
        .style("stroke-width", "0.5px");

      // bars
      let color_range = data.length <= 8 ? pals.okabe : pals.mpn65
      let color = d3.scaleOrdinal()
        .domain(data.map(d => d[groupKey]))
        .range(color_range)
      //
      svg.append("g")
        .selectAll("g")
        .data(data)
        .join("rect")
          .attr("x", d => x(0))
          .attr("y", d => y0(d[groupKey]))
          .attr("height", y0.bandwidth())
          .attr("width", d => x(d.value) - x(0))
          // .attr("fill", d => pals.okabe[0])
          .attr("fill", d => color(d[groupKey]))
          .attr("stroke", "#000000")
          .attr("stroke-width", "0.5px")

      svg.append("text")
        .attr("x", margin.left)
        .attr("y", margin.top - 5)
        .attr("font-family", "sans-serif")
        .attr("font-size", `${g_font_size}px`)
        // .attr("font-style", "italic")
        // .html(`Cells with <tspan font-style="italic">${geneSym}</tspan>`)
        .html(`Cells`)

      svg.append("g")
          .call(xAxis);

      svg.append("g")
          .call(yAxis);

      // svg.append("g")
      //     .call(legend);
    }

    function drawBars(metaInfo) {

        if (metaInfo.type != "enum" || metaInfo.valCounts.length > pals.okabe.length) {
          return
        }

        let groupKey = db.conf.clusterField
        let groups = get_groups(groupKey)

        // let subgroupKey = "case"
        // let subgroupLevels = ["Non-inflamed","Inflamed","Healthy"]
        let subgroupKey = metaInfo.name
        let subgroupLevels = metaInfo.valCounts.map(d => d[0])

        var data = []
        for (var i = 0; i < groups.length; i++) {
          var res = {
            "cluster": groups[i][0]
          }
          var counts = Array.from(d3.rollup(
            g_mydata.filter(d => d[groupKey] == groups[i][0]),
            v => v.length, d => d[subgroupKey]
          ))
          for (var j = 0; j < counts.length; j++) {
            res[counts[j][0]] = counts[j][1]
          }
          data.push(res)
        }

        let height = 500
        let width = 500
        let margin = {top: 0, right: 120, bottom: 20, left: 180}

        let y0 = d3.scaleBand()
          .domain(data.map(d => d[groupKey]))
          .rangeRound([margin.top, height - margin.bottom])
          .paddingInner(0.1)

        let y1 = d3.scaleBand()
          .domain(subgroupLevels)
          // .rangeRound([y0.bandwidth(), 0])
          .rangeRound([0, y0.bandwidth()])
          .padding(0.05)

        let x = d3.scaleLinear()
          .domain([0, d3.max(data, d => d3.max(subgroupLevels, key => d[key]))]).nice()
          .rangeRound([margin.left, width - margin.right])

        let color_range = subgroupLevels.length <= 8 ? pals.okabe : pals.mpn65
        let color = d3.scaleOrdinal()
          .range(color_range)
          // .range(["#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00"])

        let xAxis = g => g
          .attr("transform", `translate(0,${height - margin.bottom})`)
          .call(d3.axisBottom(x).tickSizeOuter(0).ticks(5, "s"))
          .call(g => g.select(".domain").remove())
          .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))
          .call(
            g => g.select(".tick:last-of-type text").clone()
              .attr("x", 15)
              .attr("text-anchor", "start")
              // .attr("font-weight", "bold")
              .attr("font-size", `${g_font_size}px`)
              .text("cells")
          )

        let yAxis = g => g
          .attr("transform", `translate(${margin.left},0)`)
          .call(d3.axisLeft(y0).ticks(null, "s"))
          .call(g => g.select(".domain").remove())
          .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))

        const legend = svg => {
          const g = svg
              .attr("transform", `translate(${width},0)`)
              .attr("text-anchor", "end")
              .attr("font-family", "sans-serif")
              .attr("font-size", `${g_font_size}px`)
            .selectAll("g")
            // .data(color.domain().slice().reverse())
            .data(color.domain().slice())
            .join("g")
              .attr("transform", (d, i) => `translate(0,${i * 20})`);
          g.append("rect")
              .attr("x", -19)
              .attr("width", 19)
              .attr("height", 19)
              .attr("fill", color);
          g.append("text")
              .attr("x", -24)
              .attr("y", 9.5)
              .attr("dy", "0.35em")
              .text(d => d);
        }

        $("#mybars").html("")
        const svg = d3.select("#mybars").append('svg')
          .attr("width", width)
          .attr("height", height)
          .style("position", "absolute")
          .style("top", '0px')
          .style("left", '0px');

        let panelBorder = svg.append("rect")
          .attr("x", margin.left)
          .attr("y", margin.top)
          .attr("height", height - margin.top - margin.bottom)
          .attr("width", width - margin.left - margin.right)
          .style("stroke", "black")
          .style("fill", "none")
          .style("stroke-width", "0.5px")

        svg.append("g")
          .selectAll("g")
          .data(data)
          .join("g")
            .attr("transform", d => `translate(0,${y0(d[groupKey])})`)
          .selectAll("rect")
          .data(d => subgroupLevels.map(key => ({key, value: d[key]})))
          .join("rect")
            .attr("x", d => x(0))
            .attr("y", d => y1(d.key))
            .attr("height", y1.bandwidth())
            .attr("width", d => x(d.value) - x(0))
            .attr("fill", d => color(d.key));

        svg.append("g")
            .call(xAxis);

        svg.append("g")
            .call(yAxis);

        svg.append("g")
            .call(legend);

        //------------------------------------------------------------------------ 
        // Vertical grouped bars
        //------------------------------------------------------------------------ 

        // let y = d3.scaleLinear()
        //   .domain([0, d3.max(data, d => d3.max(subgroupLevels, key => d[key]))]).nice()
        //   .rangeRound([height - margin.bottom, margin.top])

        // let x0 = d3.scaleBand()
        //   .domain(data.map(d => d[groupKey]))
        //   .rangeRound([margin.left, width - margin.right])
        //   .paddingInner(0.1)

        // let x1 = d3.scaleBand()
        //   .domain(subgroupLevels)
        //   .rangeRound([0, x0.bandwidth()])
        //   .padding(0.05)

        // let yAxis = g => g
        //   .attr("transform", `translate(${margin.left},0)`)
        //   .call(d3.axisLeft(y).ticks(null, "s"))
        //   .call(g => g.select(".domain").remove())
        //   .call(g => g.select(".tick:last-of-type text").clone()
        //       .attr("x", 3)
        //       .attr("text-anchor", "start")
        //       .attr("font-weight", "bold")
        //       .text(data[0]))

        // let xAxis = g => g
        //   .attr("transform", `translate(0,${height - margin.bottom})`)
        //   .call(d3.axisBottom(x0).tickSizeOuter(0))
        //   .call(g => g.select(".domain").remove())

        // let color = d3.scaleOrdinal()
        //   .range(["#98abc5", "#8a89a6", "#7b6888"])//, "#6b486b", "#a05d56", "#d0743c", "#ff8c00"])

        // $("#mybars").html("")
        // const svg = d3.select("#mybars").append('svg')
        //   .attr("width", width)
        //   .attr("height", height)
        //   .style("position", "absolute")
        //   .style("top", '0px')
        //   .style("left", '0px');

        // svg.append("g")
        //   .selectAll("g")
        //   .data(data)
        //   .join("g")
        //     .attr("transform", d => `translate(${x0(d[groupKey])},0)`)
        //   .selectAll("rect")
        //   .data(d => subgroupLevels.map(key => ({key, value: d[key]})))
        //   .join("rect")
        //     .attr("x", d => x1(d.key))
        //     .attr("y", d => y(d.value))
        //     .attr("width", x1.bandwidth())
        //     .attr("height", d => y(0) - y(d.value))
        //     .attr("fill", d => color(d.key));

        // svg.append("g")
        //     .call(xAxis);

        // svg.append("g")
        //     .call(yAxis);

        // var legend = svg => {
        //   const g = svg
        //       .attr("transform", `translate(${width},0)`)
        //       .attr("text-anchor", "end")
        //       .attr("font-family", "sans-serif")
        //       .attr("font-size", 10)
        //     .selectAll("g")
        //     .data(color.domain().slice().reverse())
        //     .join("g")
        //       .attr("transform", (d, i) => `translate(0,${i * 20})`);

        //   g.append("rect")
        //       .attr("x", -19)
        //       .attr("width", 19)
        //       .attr("height", 19)
        //       .attr("fill", color);

        //   g.append("text")
        //       .attr("x", -24)
        //       .attr("y", 9.5)
        //       .attr("dy", "0.35em")
        //       .text(d => d);
        // }

        // svg.append("g")
        //     .call(legend);

      }


    function drawGene4(data, geneSym, gene_groupby) {
      console.log(`drawGene4(data, geneSym = ${geneSym}, gene_groupby = ${gene_groupby})`)
      console.log({data})
      //
      let canvas_scale = 2
      let radius = g_size.hex_radius
      let legend_width = 90
      let plot_width = 360 * canvas_scale
      let plot_height = 360 * canvas_scale
      //
      let legend_margin = {
        top: 70, right: 0, bottom: 0, left: 5
      }
      let margin = {
        top: 55, right: 10, bottom: 10, left: 10
      }
      //
      var aggKey = g_donorField
      var metaInfo = null
      var subgroupKey = null
      var subgroupLevels = []
      var subgroupCounts = {}
      if (gene_groupby !== "none") {
        metaInfo = db.findMetaInfo(gene_groupby)
        if (metaInfo.type === "enum" || metaInfo.valCounts.length < 5) {
          subgroupKey = metaInfo.name
          subgroupLevels = metaInfo.valCounts.map(d => d[0]).slice().sort()
          subgroupCounts = Object.fromEntries(d3.rollup(
            data, v => [...new Set(v.map(d => d[aggKey]))].length, d => d[subgroupKey]
          ))
        }
      }
      //
      d3.select("#gene-loading-spinner").remove()
      d3.select("#embedding-gene")
        .html("")
        .append("div")
        .attr("class", "row")
        .attr("id", "gene-row")
      //
      const plot_total = subgroupKey ? subgroupLevels.length : 1
      if (plot_total === 1) {
        margin.top = 35
      }
      if (plot_total > 1) {
        radius = radius * 3
      }
      var plot_col_width = 12
      let font_size = g_font_size * 1.2
      if (plot_total > 1) {
        plot_col_width = 6
        font_size = 24
      }
      //
      var x = d3.scaleLinear()
        .domain(d3.extent(data, d => d.x))
        .range([margin.left * canvas_scale, plot_width - margin.right * canvas_scale])

      var y = d3.scaleLinear()
        .domain(d3.extent(data, d => d.y))
        .range([-30 + plot_height - margin.bottom * canvas_scale, margin.top * canvas_scale])

      var hexbin = d3.hexbin()
        .x(d => x(d.x))
        .y(d => y(d.y))
        .radius(radius * (plot_width - legend_width * canvas_scale) / (plot_height))
        .extent([
          [margin.left, margin.top],
          [plot_width - margin.right, plot_height - margin.bottom]
        ])
      var bins = hexbin(data)
      // Sequential scale (works)
      var bin_max = d3.max(bins, bin => d3.mean(bin, d => d.gene))
      // var color = d3.scaleSequential(d3.interpolateBuPu)
      //   .domain([0, bin_max])
      var color = d3.scaleLinear()
        .domain(linspace(0, bin_max, 100))
        .range(pals.batlow)
      var color_range = linspace(0, bin_max, 10).reverse()
      var hex = new Path2D(hexbin.hexagon())
      //
      for (var plot_i = 0; plot_i < plot_total; plot_i++) {
        //
        const div_id =`embedding-gene-plot-${plot_i}`
        const container = d3.select("#gene-row")
          .append("div")
          .attr("class", `col-${plot_col_width}`)
          .attr("id", div_id)
        //
        const canvas_id = `mycanvas-gene-plot-${plot_i}`
        d3.select(canvas_id).remove()
        const canvas = container.append('canvas')
          .attr("id", canvas_id)
          .style("width", "100%")
          .node()
        canvas.width = plot_width + legend_width * canvas_scale
        canvas.height = plot_height

        const context = canvas.getContext('2d')
        context.fillStyle = "#fff"

        const svg = container.append('svg')
          .attr("viewBox",
            `0 0 ${canvas.width / canvas_scale} ${1 + canvas.height / canvas_scale}`
          )
          .style("position", "absolute")
          .style("left", "10px")
          .style("top", "0px")
          .style("height", "100%")
          // .style("top", $(`#${div_id}`).css("padding-top"))
          // .style("left", $(`#${div_id}`).css("padding-left"))

        if (plot_i === 0 && plot_total > 1) {
          svg.append("text")
            .attr("x", margin.left)
            .attr("y", margin.top - 35)
            .attr("font-family", "sans-serif")
            .attr("font-size", `${font_size * 1.1}px`)
            .attr("font-style", "italic")
            .text(geneSym)
        } else if (plot_i === 0 && plot_total === 1) {
          svg.append("text")
            .attr("x", margin.left)
            .attr("y", margin.top - 15)
            .attr("font-family", "sans-serif")
            .attr("font-size", `${font_size * 1.4}px`)
            .attr("font-style", "italic")
            .text(geneSym)
        }
        
        // Count the cells per subgroup
        if (plot_total === 1) {
          var n_cells = data.filter(d => d.gene > 0).length
          var total_cells = data.length
        } else {
          var n_cells = data
            .filter(
              d => d.gene > 0 && d[subgroupKey] == subgroupLevels[plot_i]
            ).length
          var total_cells = data
            .filter(d => d[subgroupKey] == subgroupLevels[plot_i])
            .length
        }

        if (plot_total > 1) {
          // var subgroup_n = subgroupCounts[subgroupLevels[plot_i]]
          svg.append("text")
            .attr("x", (plot_width / canvas_scale) / 2)
            .attr("y", margin.top - 15)
            .attr("text-anchor", "middle")
            .attr("font-family", "sans-serif")
            .attr("font-size", `${font_size}px`)
            .text(`${subgroupLevels[plot_i]} (n = ${d3.format(",")(total_cells)})`)
        }

        svg.append("text")
          .attr("x", margin.left / 2)
          .attr("y", plot_height / canvas_scale - margin.bottom / 2)
          .attr("font-family", "sans-serif")
          .attr("font-size", `${font_size}px`)
          .text(`${d3.format(",")(n_cells)} (${d3.format(".1%")(n_cells / total_cells)}) cells`)

        let panelBorder = svg.append("rect")
          .attr("x", margin.left - 9)
          .attr("y", margin.top - 10)
          .attr("width", plot_width / canvas_scale)
          .attr("height", plot_height / canvas_scale - margin.top + margin.bottom)
          .style("stroke", "black")
          .style("fill", "none")
          .style("stroke-width", "0.5px");

        bins.forEach(function(bin) {
          var color_val = 0
          var my_bin = bin.filter(d => d[subgroupKey] == subgroupLevels[plot_i])
          if (my_bin.length) {
            if (plot_total > 1) {
              color_val = d3.mean(
                my_bin,
                d => d.gene
              )
            } else {
              color_val = d3.mean(bin, d => d.gene)
            }
            context.translate(bin.x, bin.y);
            //context.fillStyle = color(bin.length);
            context.fillStyle = color(color_val)
            context.fill(hex)
            context.strokeStyle = context.fillStyle
            context.lineWidth = 0
            context.stroke(hex)
            context.setTransform(1, 0, 0, 1, 0, 0)
          }
        })

        if (plot_total === 1 || plot_i === plot_total - 1) {
          svg.append("text")
            .attr("x", plot_width / canvas_scale + legend_margin.left)
            .attr("y", legend_margin.top - 20)
            .attr("font-family", "sans-serif")
            .attr("font-size", `${font_size}px`)
            .attr("font-weight", "bold")
            .text("log2CPM")
          var legend = svg.selectAll("g.legend_colorbar")
            .data(color_range)
            .enter()
            .append("g")
            .attr("class","legend_colorbar");
          let legend_rect_height = 25
          legend
            .append('rect')
            .attr("x", plot_width / canvas_scale + legend_margin.left)
            .attr("y", function(d, i) {
               return legend_margin.top + i * legend_rect_height;
            })
           .attr("width", 15)
           .attr("height", legend_rect_height)
           .style("fill", function(d){return color(d);});
          legend
            .append('text')
            .attr("x", plot_width / canvas_scale + legend_margin.left + 20)
            .attr("y", function(d, i) {
             return legend_margin.top + i * legend_rect_height;
            })
            .attr("alignment-baseline", "middle")
            .style("font-size", `${font_size}px`)
            .text(function(d, i) {
              var format = d3.format(".1f");
              return `${format(+color_range[i])}`
            })
        }
        if (plot_total == 2 && plot_i == 1) {
          const blank_container = d3.select("#gene-row")
            .append("div")
            .attr("class", `col-${plot_col_width}`)
            // .style("position", "relative")
            // .style("display", "inline-block")
          const blank_canvas = blank_container.append('canvas')
            .attr("id", "gene-blank")
            .style("width", "100%")
            .node()
          blank_canvas.width = plot_width + legend_width * canvas_scale
          blank_canvas.height = plot_height
        }
      } // for plot_i
      // layout hack
      if (plot_total > 1) {
        d3.select("#embedding-gene").attr("class", "col-12 mb-3")
      } else {
        d3.select("#embedding-gene").attr("class", "col-12")
      }
    }

    function draw_meta_boxplot(fieldName, meta_groupby) {
      if (meta_groupby === "none") {
        return draw_meta_boxplot_none()
      }
      //
      let data = g_mydata
      const groupKey = db.conf.clusterField
      const fillKey = meta_groupby // TBStatus, case/control, etc.
      const aggKey = g_donorField
      // const agg_counts = Object.fromEntries(
      //   db.conf.metaFields.filter(d => d.name == aggKey)[0].valCounts
      // )
      const agg_counts = Object.fromEntries(get_groups(aggKey))
      let metaInfo = db.findMetaInfo(fillKey)
      let subgroupLevels = metaInfo.valCounts.map(d => d[0]).slice().sort()
      let subgroupCounts = Object.fromEntries(d3.rollup(
        data, v => [...new Set(v.map(d => d[aggKey]))].length, d => d[fillKey]
      ))
      var make_bin = function(d) {
        // percent of each donor's cells in each cluster
        const values = Array.from(d3.rollup(
          d, v => v.length, d => d[aggKey]
        )).map(d => d[1] / agg_counts[d[0]])
          .sort((a, b) => a - b)
        const min = values[0]
        const max = values[values.length - 1]
        const q1 = d3.quantile(values, 0.25)
        const q2 = d3.quantile(values, 0.5)
        const q3 = d3.quantile(values, 0.75)
        const iqr = q3 - q1
        let bin = {
          min: min,
          max: max,
          q1: q1,
          q2: q2,
          q3: q3,
          iqr: q3 - q1,
          r0: Math.max(min, q1 - iqr * 1.5),
          r1: Math.min(max, q3 + iqr * 1.5)
        }
        bin.outliers = values.filter(v => v < bin.r0 || v > bin.r1)
        return bin 
      }
      //
      let bins = []
      Array.from(d3.group(data, d => d[groupKey])).map(d => {
        const m = d3.group(d[1], v => v[fillKey])
        for (const [k, v] of m.entries()) {
          var bin = {}
          bin[groupKey] = d[0]
          bin[fillKey] = k
          bin["bin"] = make_bin(v)
          bins.push(bin)
        }
      })
      bins = bins //.sort((a,b) => a[groupKey].localeCompare(b[groupKey]))
      //
      let groups = get_groups(groupKey)
      let group_names = groups.map(d => d[0])
      var el = $("#meta-boxplot")
      el.html("")
      const longest_subgroup_key = d3.greatest(
        subgroupLevels, (a, b) => d3.ascending(a.length, b.length)
      )
      const longest_group_key = d3.greatest(
        group_names, (a, b) => d3.ascending(a.length, b.length)
      )
      const margin = {
        top: 15,
        right: getTextWidth(longest_subgroup_key, "15px arial") + 80,
        bottom: 20,
        left: getTextWidth(longest_group_key, "15px arial") + 9
      }
      const width = 500
      const height = Math.max(24 * groups.length * (0.7 * subgroupLevels.length), 300)
      if (height > el.outerHeight()) {
        d3.select("#meta-row").style("height", `${height + margin.top + margin.bottom}px`)
      }
      //
      const svg = d3.select("#meta-boxplot").append('svg')
        // .attr("width", width)
        .attr("height", "100%")
        .attr("viewBox", `0 0 ${width} ${height}`)
        // .style("position", "absolute")
        // .style("top", '0px')
        // .style("left", '0px');
      //
      var y0 = d3.scaleBand()
        .domain(group_names)
        // .domain(bins.map(d => d[groupKey]))
        // .rangeRound([margin.top, height - margin.bottom])
        .range([margin.top + 2, height - margin.bottom - 2])
        .paddingInner(0.1)
        .paddingOuter(0.0)
      //
      var y1 = d3.scaleBand()
        .domain(subgroupLevels)
        .rangeRound([0, y0.bandwidth()])
        .padding(0.05)
      //
      var x = d3.scaleLog()
      // var x = d3.scaleLinear()
        .domain([d3.min(bins, d => d.bin.min), d3.max(bins, d => d.bin.max)]).nice()
        .rangeRound([margin.left + 5, width - margin.right - 5])
      //
      var xAxis = g => g
        .attr("transform", `translate(0,${height - margin.bottom})`)
        // .call(d3.axisBottom(x).tickSizeOuter(0).ticks(5, "~d"))
        .call(d3.axisBottom(x).tickSizeOuter(0).ticks(5, "~p"))
        .call(g => g.select(".domain").remove())
        .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))
        // .call(d3.axisBottom(x).ticks(10, [5, "~g"]))
        // .call(d3.axisBottom(x).ticks(10, [5, "~d"]))
        // .call(g => g.select(".tick:last-of-type text").clone()
        //     .attr("x", 15)
        //     .attr("y", -5)
        //     .attr("text-anchor", "start")
        //     .attr("font-weight", "bold")
        //     .attr("font-size", 14)
        //     .text("Percent"))
      //
      var yAxis = g => g
        .attr("transform", `translate(${margin.left},0)`)
        .call(d3.axisLeft(y0).ticks(null, "s"))
        .call(g => g.select(".domain").remove())
        .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))
        .call(g => g.selectAll(".tick text") .text(d => d))
      //
      let color_range = subgroupLevels.length <= 8 ? pals.okabe : pals.mpn65
      if (meta_groupby === 'case') {
        color_range = pals.okabe.slice(0,2).reverse()
      }
      let color = d3.scaleOrdinal()
        .domain(subgroupLevels)
        .range(color_range)
      // background stripes
      const bg = svg.append("g")
        .selectAll("g")
        .data(groups)
        .join("g")
      bg.append("rect")
        .attr("fill", (d, i) => i % 2 == 0 ? "#ffffff00" : "#eee")
        .attr("x", margin.left)
        .attr("y", d => y0(d[0]) - y1.bandwidth() / 8)
        .attr("height", y0.bandwidth() + y1.bandwidth() / 4)
        .attr("width", width - margin.left - margin.right)
      // vertical lines
      const vLines = svg.append("g")
        .selectAll("g")
        .data([0.001, 0.01, 0.1, 1.0])
        .join("g")
      vLines.append("path")
        // .attr("stroke", "#eee")
        .attr("stroke", "#fff")
        .attr("stroke-width", "2px")
        .attr("d", d => `M${x(d)},${margin.top} V${height - margin.bottom}`);
      // let vLines = svg.append("g")
      //   .style("stroke", "#000")
      //   .style("stroke-width", "0.5px")
      //   .append("path")
      //   .attr("d", `M${x(0.1)},${margin.top} V${height - margin.bottom}`)
      //
      const g = svg.append("g")
        .selectAll("g")
        .data(bins)
        .join("g")
      // range lines
      g.append("path")
          .attr("transform",
            d => `translate(0,${y0(d[groupKey]) + y1(d[fillKey]) + y1.bandwidth()})`
          )
          .attr("stroke", "currentColor")
          .attr("stroke-width", "0.5px")
          .attr("d", d => `
            M${x(d.bin.r0)},${-y1.bandwidth() / 2}
            H${x(d.bin.r1)}
          `)
      g.append("g")
          .attr("transform",
            d => `translate(0,${y0(d[groupKey]) + y1(d[fillKey]) + y1.bandwidth() / 2})`
          )
        .selectAll("circle")
        .data(d => d.bin.outliers)
        .join("circle")
          .attr("cx", d => x(d))
          .attr("r", 1)
      // boxes
      g.append("path")
          .attr("fill", "#333")
          .attr("transform",
            d => `translate(0,${y0(d[groupKey]) + y1(d[fillKey]) + y1.bandwidth()})`
          )
          .attr("fill", d => color(d[fillKey]))
          .attr("stroke", "#000000")
          .attr("stroke-width", "0.5px")
          .attr("d", d => `
            M${x(d.bin.q3)},0
            V${-y1.bandwidth()}
            H${x(d.bin.q1)}
            V${0}
            Z
          `);
      // median lines
      g.append("path")
          .attr("stroke", "currentColor")
          .attr("stroke-width", 2)
          .attr("d", d => `
            M${x(d.bin.q2)},${y0(d[groupKey]) + y1(d[fillKey]) + y1.bandwidth()}
            V${y0(d[groupKey]) + y1(d[fillKey])}
          `)
      //
      var legend = svg => {
        const g = svg
            .attr("transform", `translate(${width},${margin.top})`)
            .attr("text-anchor", "end")
            .attr("font-family", "sans-serif")
            .attr("font-size", `${g_font_size}px`)
          .selectAll("g")
          .data(color.domain())
          .join("g")
            .attr("transform", (d, i) => `translate(0,${i * 20})`);
        g.append("rect")
            .attr("x", -19)
            .attr("width", 19)
            .attr("height", 19)
            .attr("fill", color);
        g.append("text")
            .attr("x", -24)
            .attr("y", 9.5)
            .attr("dy", "0.35em")
            .text(d => `${d} (n = ${subgroupCounts[d]})`);
      }
      //
      let panelBorder = svg.append("rect")
        .attr("x", margin.left)
        .attr("y", margin.top)
        .attr("height", height - margin.top - margin.bottom)
        .attr("width", width - margin.left - margin.right)
        .style("stroke", "black")
        .style("fill", "none")
        .style("stroke-width", "0.5px");
      //
      svg.append("text")
        .attr("x", margin.left)
        .attr("y", margin.top - 5)
        .attr("font-family", "sans-serif")
        .attr("font-size", `${g_font_size}px`)
        // .attr("font-style", "italic")
        .html(`Percent of each donor's cells`)
      // draw
      svg.append("g")
          .call(xAxis);
      svg.append("g")
          .call(yAxis);
      svg.append("g")
          .call(legend);
    }

    function get_groups(groupKey) {
      return db.conf.metaFields
        .filter(d => d.name == groupKey)[0]
        .valCounts.slice()
        .sort((a, b) => a[0].localeCompare(b[0], navigator.languages[0] || navigator.language, {numeric: true, ignorePunctuation: true}))
    }

    function draw_meta_heatmap(groupKey, subgroupKey) {
      $("#meta-heatmap").html("")
      //
      let metaInfo = db.findMetaInfo(subgroupKey)
      let subgroupLevels = metaInfo.valCounts.map(d => d[0]).slice().sort()
      //
      const aggKey = g_donorField
      let groups = get_groups(groupKey)
      let group_names = groups.map(d => d[0])
      let data = Array.from(d3.rollup(
        g_mydata,
        v => Object.fromEntries(Array.from(
          d3.rollup(v, v => v.length, d => d[groupKey])
        )),
        d => [ d[aggKey], d[subgroupKey] ].join('\t')
        // d => d[aggKey]
      )).map(d => ({
        name: d[0],
        cells: d[1]
      }))
      var color_domain = []
      for (var i = 0; i < data.length; i++) {
        var cells = []
        for (var j = 0; j < groups.length; j++) {
          var group_name = groups[j][0]
          if (group_name in data[i].cells) {
            var val = data[i].cells[group_name]
            cells.push(val)
            color_domain.push(val)
          } else {
            cells.push(0)
            color_domain.push(0)
          }
        }
        data[i].cells = cells
      }
      var hc = hcluster()
        .distance('euclidean') // support for 'euclidean' and 'angular'
        .linkage('avg')        // support for 'avg', 'max' and 'min'
        .posKey('cells')    // 'position' by default
        .data(data)         // as an array of objects w/ array values for 'position'
      data = hc.orderedNodes()
      //
      var el = $("#meta-heatmap")
      el.html("")
      const legend_margin = {top: 50, right: 10, bottom: 10, left: 10}
      // const width = Math.min(1200, 8 * data.length)
      const width = el.outerWidth()
      const height = Math.max(20 * (groups.length + 2), el.outerHeight())
      const longest_group_key = d3.greatest(
        group_names, (a, b) => d3.ascending(a.length, b.length)
      )
      const margin = {
        top: 50, right: 160, bottom: 50,
        left: getTextWidth(longest_group_key, "15px arial") + 9
      }
      d3.select("#meta-row-2").style("height", `${height + margin.top + margin.bottom}px`)
      const svg = d3.select("#meta-heatmap").append('svg')
        .attr("width", width)
        .attr("height", height)
        .style("position", "absolute")
        .style("top", '0px')
        .style("left", '0px');
      // cluster
      var y = d3.scaleBand()
        .domain(groups.map(d => d[0]))
        .range([margin.top + 2, height - margin.bottom - 2])
        .paddingInner(0.05)
        // .paddingOuter(0.0)
      // donor
      var x = d3.scaleBand()
        .domain(data.map(d => d.name))
        .range([margin.left + 2, width - margin.right - 2])
      //
      var xAxis = g => g
        .attr("transform", `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(null, "s"))
        .call(g => g.select(".domain").remove())
        .call(
          g => g.selectAll(".tick text")
            .attr("font-size", `${g_font_size}px`)
            .attr("transform", "rotate(-45) translate(-7,-8)")
            .attr("text-anchor", "end")
            .text(d => d.split('\t')[0])
        )
      //
      var yAxis = g => g
        .attr("transform", `translate(${margin.left},0)`)
        .call(d3.axisLeft(y).ticks(null, "s"))
        .call(g => g.select(".domain").remove())
        .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))
      //
      let color_range = linspace(0, 1, 11)
      var color = d3.scaleQuantile()
        .domain(color_domain)
        .range(color_range.map(t => d3.interpolatePuBuGn(t)))
      //
      let color_range_subgroup = subgroupLevels.length <= 8 ? pals.okabe : pals.mpn65
      let color_subgroup = d3.scaleOrdinal()
        .domain(subgroupLevels)
        .range(color_range_subgroup)
      const g = svg.append("g")
        .selectAll("g")
        .data(data)
        .join("g")
          .attr("transform", (d, i) => `translate(${x(d.name)},0)`)
        .selectAll("rect")
        .data(d => d.cells)
        .join("rect")
          .attr("x", 0)
          .attr("y", (d, i) => y(group_names[i]))
          .attr("width", x.bandwidth())
          .attr("height", y.bandwidth())
          .attr("fill", d => color(d))
      // subgroup bar
      svg.append("g")
        .selectAll("rect")
        .data(data)
        .join("rect")
          .attr("transform", (d, i) => `translate(${x(d.name)},0)`)
          .attr("x", 0)
          .attr("y", margin.top - y.bandwidth() + 1)
          .attr("width", x.bandwidth())
          .attr("height", y.bandwidth() - 3)
          .attr("fill", d => color_subgroup(d.name.split('\t')[1]))
      //
      svg.append("rect")
        .attr("x", margin.left)
        .attr("y", margin.top - y.bandwidth() - 1)
        .attr("height", y.bandwidth() + 1)
        .attr("width", width - margin.left - margin.right)
        .style("stroke", "black")
        .style("fill", "none")
        .style("stroke-width", "1px");
      //
      let panelBorder = svg.append("rect")
        .attr("x", margin.left)
        .attr("y", margin.top)
        .attr("height", height - margin.top - margin.bottom)
        .attr("width", width - margin.left - margin.right)
        .style("stroke", "black")
        .style("fill", "none")
        .style("stroke-width", "0.5px");
      //
      svg.append("text")
        .attr("x", width - margin.right + legend_margin.left)
        .attr("y", margin.top - y.bandwidth() / 2)
        .attr("font-family", "sans-serif")
        .attr("font-size", `${g_font_size}px`)
        .attr("font-weight", "bold")
        // .attr("alignment-baseline", "after-edge")
        .attr("alignment-baseline", "middle")
        .html(subgroupKey)
      //
      let n_agg = [...new Set(data.map(d => d.name.split('\t')[0]))].length
      svg.append("text")
        .attr("x", margin.left)
        .attr("y", margin.top - y.bandwidth() - 5)
        .attr("font-family", "sans-serif")
        .attr("font-size", `${g_font_size}px`)
        // .attr("font-style", "italic")
        .html(`Number of cells in each ${groupKey}, per ${aggKey} (n = ${n_agg})`)
      // draw
      if (data.length < 50) {
        svg.append("g")
            .call(xAxis)
      }
      svg.append("g")
          .call(yAxis)
      // legend
      // svg.append("text")
      //   .attr("x", width - margin.right + legend_margin.left)
      //   .attr("y", legend_margin.top)
      //   .attr("font-family", "sans-serif")
      //   .attr("font-size", `${g_font_size}px`)
      //   .attr("font-weight", "bold")
      //   .text("Cells")
      var legend = svg.selectAll("g.legend_colorbar")
        .data(color.range().reverse())
        .enter()
        .append("g")
        .attr("class","legend_colorbar")
        .attr("transform",
          `translate(${width - margin.right + legend_margin.left},${legend_margin.top})`
        )
      const legend_rect = {width: x.bandwidth() / 1.5, height: y.bandwidth()}
      legend
        .append('rect')
        .attr("y", (d, i) => (i + 1) * legend_rect.height)
        .attr("width", legend_rect.width)
        .attr("height", legend_rect.height)
        // .style("stroke", "black")
        // .style("stroke-width", 0.1)
        .style("fill", d => d)
      legend
        .append('text')
        .attr("x", legend_rect.width + 2)
        .attr("y", (d, i) => (i + 1) * legend_rect.height)
        .attr("alignment-baseline", "middle")
        .style("font-size", `${g_font_size}px`)
        .text(function(d, i) {
          var extent = color.invertExtent(d)
          // var format = d3.format(".2")
          return `${d3.format("d")(+extent[1])}`
        })
      legend
        .append("text")
        .attr("x", legend_rect.width + 2)
        .attr("y", (d, i) => (i + 2) * legend_rect.height)
        .attr("alignment-baseline", "middle")
        .style("font-size",`${g_font_size}px`)
        .text(function(d, i) {
          if (i == color_range.length - 1) {
            var extent = color.invertExtent(d)
            return `${d3.format("d")(+extent[0])}`
          }
        })
    }

    function draw_gene_boxplot(geneSym, gene_groupby) {
      if (gene_groupby === "none") {
        return draw_gene_boxplot_none(geneSym)
      }
      //
      let data = g_mydata
      const groupKey = db.conf.clusterField
      const fillKey = gene_groupby
      const aggKey = g_donorField
      let metaInfo = db.findMetaInfo(fillKey)
      let subgroupLevels = metaInfo.valCounts.map(d => d[0]).slice().sort()
      let subgroupCounts = Object.fromEntries(d3.rollup(
        data, v => [...new Set(v.map(d => d[aggKey]))].length, d => d[fillKey]
      ))
      const longest_group_key = d3.greatest(
        Object.keys(subgroupCounts), (a, b) => d3.ascending(a.length, b.length)
      )
      var make_bin = function(d) {
        // d.sort((a, b) => a.gene - b.gene)
        // const values = d.map(d => d.gene).filter(x => x > 0).sort()
        // mean by donor
        var values = Object.values(Object.fromEntries(d3.rollup(
          d, v => d3.mean(v.map(v => v.gene)), d => d[aggKey]
        ))).sort((a,b) => a - b)
        if (values.length > 2) {
          const min = values[0]
          const max = values[values.length - 1]
          const q1 = d3.quantile(values, 0.25)
          const q2 = d3.quantile(values, 0.5)
          const q3 = d3.quantile(values, 0.75)
          const iqr = q3 - q1
          let bin = {
            min: min,
            max: max,
            q1: q1,
            q2: q2,
            q3: q3,
            iqr: q3 - q1,
            r0: Math.max(min, q1 - iqr * 1.5),
            r1: Math.min(max, q3 + iqr * 1.5)
          }
          bin["outliers"] = values.filter(v => v < bin.r0 || v > bin.r1)
          return bin 
        }
        return {
          min: 0,
          max: 0,
          q1: 0,
          q2: 0,
          q3: 0,
          iqr: 0,
          r0: 0,
          r1: 0,
          outliers: []
        }
      }
      //
      let bins = []
      Array.from(d3.group(data, d => d[groupKey])).map(d => {
        const m = d3.group(d[1], v => v[fillKey])
        for (const [k, v] of m.entries()) {
          var bin = {}
          bin[groupKey] = d[0]
          bin[fillKey] = k
          bin["bin"] = make_bin(v)
          bins.push(bin)
        }
      })
      // bins = bins.sort((a,b) => a[groupKey].localeCompare(b[groupKey]))
      //
      let groups = get_groups(groupKey)
      let group_names = groups.map(d => d[0])
      //
      var el = $("#gene-boxplot")
      el.html("")
      const width = g_size.gene_boxplot.width
      const height = g_size.gene_boxplot.height(groups.length, subgroupLevels.length)
      const margin = g_size.gene_boxplot.margin
      // const margin = {
      //   top: 15,
      //   right: 5,
      //   // right: getTextWidth(longest_group_key, "14px arial") + 80,
      //   // right: 165,
      //   bottom: 40,
      //   left: 5,
      //   xoffset: g_size.gene_bars.width
      // }
      //
      // const svg = d3.select("#gene-boxplot").append('svg')
      //   .attr("viewBox", `0 0 ${width} ${height}`)
      //   .style("position", "absolute")
      //   .style("top", '0px')
      //   .style("left", '0px')
      const svg = d3.select("#gene-bars-svg").append("g")
        .attr("transform", `translate(${margin.xoffset},0)`)
      //
      var y0 = d3.scaleBand()
        .domain(group_names)
        // .domain(bins.map(d => d[groupKey]))
        // .rangeRound([margin.top, height - margin.bottom])
        .range([margin.top + 2, height - margin.bottom - 2])
        .paddingInner(0.1)
        .paddingOuter(0.0)
      //
      // console.log("boxplot")
      // console.log(subgroupLevels)
      var y1 = d3.scaleBand()
        .domain(subgroupLevels)
        .rangeRound([0, y0.bandwidth()])
        .padding(0.05)
      //
      // var x = d3.scaleLog()
      var x = d3.scaleLinear()
        // .domain([d3.min(data, d => d.gene), d3.max(data, d => d.gene)]).nice()
        .domain([d3.min(bins, d => d.bin.min), d3.max(bins, d => d.bin.max)])
        .rangeRound([margin.left + 5, width - margin.right - 5])
      //
      var xAxis = g => g
        .attr("transform", `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).tickSizeOuter(0).ticks(5, "~g"))
        .call(g => g.select(".domain").remove())
        .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))
      // x-axis label right
      svg.append("text")
        // .attr("x", width - margin.right + 5)
        // .attr("y", height - margin.bottom)
        .attr("x", (width / 2) - 30)
        .attr("y", height - margin.bottom / 5)
        .attr("alignment-baseline", "middle")
        .attr("font-family", "sans-serif")
        .attr("font-size", `${g_font_size}px`)
        .attr("font-weight", "bold")
        .text("Log2CPM")
      //
      var yAxis = g => g
        .attr("transform", `translate(${margin.left},0)`)
        .call(d3.axisLeft(y0).ticks(null, "s"))
        .call(g => g.select(".domain").remove())
        .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))
      //
      let color_range = subgroupLevels.length <= 8 ? pals.okabe : pals.mpn65
      if (gene_groupby === 'case') {
        color_range = pals.okabe.slice(0,2).reverse()
      }
      let color = d3.scaleOrdinal()
        .domain(subgroupLevels)
        .range(color_range)
      // background stripes
      const bg = svg.append("g")
        .selectAll("g")
        .data(groups)
        .join("g")
      bg.append("rect")
        .attr("fill", (d, i) => i % 2 == 0 ? "#ffffff00" : "#eee")
        .attr("x", margin.left)
        .attr("width", width - margin.left - margin.right)
        .attr("y", d => y0(d[0]))
        .attr("height", y0.bandwidth())
      // 
      const g = svg.append("g")
        .selectAll("g")
        .data(bins)
        .join("g")
      // range lines
      g.append("path")
          .attr("transform",
            d => `translate(0,${y0(d[groupKey]) + y1(d[fillKey]) + y1.bandwidth()})`
          )
          .attr("stroke", "currentColor")
          .attr("stroke-width", "0.5px")
          .attr("d", d => `
            M${x(d.bin.r0)},${-y1.bandwidth() / 2}
            H${x(d.bin.r1)}
          `);
      g.append("g")
          .attr("transform",
            d => `translate(0,${y0(d[groupKey]) + y1(d[fillKey]) + y1.bandwidth() / 2})`
          )
        .selectAll("circle")
        .data(d => d.bin.outliers)
        .join("circle")
          .attr("cx", d => x(d))
          .attr("r", 1)
      // boxes
      g.append("path")
          .attr("fill", "#333")
          .attr("transform",
            d => `translate(0,${y0(d[groupKey]) + y1(d[fillKey]) + y1.bandwidth()})`
          )
          .attr("fill", d => color(d[fillKey]))
          .attr("stroke", "#000000")
          .attr("stroke-width", "0.5px")
          .attr("d", d => `
            M${x(d.bin.q3)},0
            V${-y1.bandwidth()}
            H${x(d.bin.q1)}
            V${0}
            Z
          `);
      // median lines
      g.append("path")
          .attr("stroke", "currentColor")
          .attr("stroke-width", 2)
          .attr("d", d => `
            M${x(d.bin.q2)},${y0(d[groupKey]) + y1(d[fillKey]) + y1.bandwidth()}
            V${y0(d[groupKey]) + y1(d[fillKey])}
          `);
      //
      let panelBorder = svg.append("rect")
        .attr("x", margin.left)
        .attr("y", margin.top)
        .attr("height", height - margin.top - margin.bottom)
        .attr("width", width - margin.left - margin.right)
        .style("stroke", "black")
        .style("fill", "none")
        .style("stroke-width", "0.5px");
      //
      svg.append("text")
        .attr("x", margin.left)
        .attr("y", margin.top - 5)
        .attr("font-family", "sans-serif")
        .attr("font-size", `${g_font_size}px`)
        // .attr("font-style", "italic")
        .html(`Mean <tspan font-style="italic">${geneSym}</tspan> by ${aggKey}`)
      // draw
      svg.append("g")
          .call(xAxis);
      // svg.append("g")
      //     .call(yAxis);
      // draw legend separately
      // svg.append("g")
      //     .call(legend);
      $("#gene-legend").html("")
      //
      var legend_text_width = getTextWidth(longest_group_key + ' (n = 100)', `${g_font_size}px arial`)
      var legend = svg => {
        const g = svg
            // .attr("transform", `translate(${width},${margin.top})`)
            .attr("transform", `translate(${legend_text_width},0)`)
            .attr("text-anchor", "end")
            .attr("font-family", "sans-serif")
            .attr("font-size", `${g_font_size}px`)
          .selectAll("g")
          .data(color.domain())
          .join("g")
            .attr("transform", (d, i) => `translate(0,${i * 20})`)
        g.append("rect")
            .attr("x", 2)
            .attr("width", 16)
            .attr("height", 16)
            .attr("fill", color)
        g.append("text")
            .attr("x", 0)
            .attr("y", 8)
            .attr("dy", "0.35em")
            .text(d => `${d} (n = ${subgroupCounts[d]})`);
            // .text(d => d);
      }
      const svg_legend = d3.select("#gene-legend").append('svg')
        .attr("viewBox", `0 0 ${legend_text_width + 20} ${(1 + color.domain().length) * 20}`)
        // .attr("width", width)
        .attr("height", (1 + color.domain().length) * 24)
        .style("float", "right")
        // .style("position", "absolute")
        // .style("top", '0px')
        // .style("left", '25px')
        // .style("left", '-25px')
       svg_legend.append("g").call(legend)
       // d3.select("#gene-legend").style("min-height", `${height}px;`)
    }

    function draw_gene_boxplot_none(geneSym) {
      let data = g_mydata
      //
      const groupKey = db.conf.clusterField
      if (!groupKey in data[0]) {
        console.log(`Object.keys(data[0]) is ${Object.keys(data[0])}`)
        return
      }
      // console.log(`draw_gene_boxplot_none('${geneSym}') groupKey = ${groupKey}`)
      const aggKey = g_donorField
      var make_bin = function(d) {
        // d.sort((a, b) => a.gene - b.gene)
        // const values = d.map(d => d.gene)
        //   .filter(x => x > 0)
        //   .sort((a, b) => a - b)
        // mean by donor
        var values = Object.values(Object.fromEntries(d3.rollup(
          d, v => d3.mean(v.map(v => v.gene)), d => d[aggKey]
        ))).sort((a,b) => a - b)
        if (values.length > 2) {
          const min = values[0]
          const max = values[values.length - 1]
          const q1 = d3.quantile(values, 0.25)
          const q2 = d3.quantile(values, 0.5)
          const q3 = d3.quantile(values, 0.75)
          const iqr = q3 - q1
          let bin = {
            min: min,
            max: max,
            q1: q1,
            q2: q2,
            q3: q3,
            iqr: q3 - q1,
            r0: Math.max(min, q1 - iqr * 1.5),
            r1: Math.min(max, q3 + iqr * 1.5)
          }
          bin["outliers"] = values.filter(v => v < bin.r0 || v > bin.r1)
          return bin 
        }
        return {
          min: 0,
          max: 0,
          q1: 0,
          q2: 0,
          q3: 0,
          iqr: 0,
          r0: 0,
          r1: 0,
          outliers: []
        }
      }
      //
      let bins = []
      Array.from(d3.group(data, d => d[groupKey])).map(d => {
        var bin = {}
        bin[groupKey] = d[0]
        bin["bin"] = make_bin(d[1])
        bins.push(bin)
      })
      bins = bins.sort((a, b) => a[groupKey].localeCompare(b[groupKey], navigator.languages[0] || navigator.language, {numeric: true, ignorePunctuation: true}))
      //
      let groups = get_groups(groupKey)
      let group_names = groups.map(d => d[0])
      //
      var el = $("#gene-boxplot")
      el.html("")
      $("#gene-legend").html("")
      const width = g_size.gene_boxplot.width
      const height = g_size.gene_boxplot.height(groups.length, 1)
      const margin = g_size.gene_boxplot.margin
      //
      // const svg = d3.select("#gene-boxplot").append('svg')
      //   .attr("viewBox", `0 0 ${width} ${height}`)
      //   .style("position", "absolute")
      //   .style("top", '0px')
      //   .style("left", '0px');
      const svg = d3.select("#gene-bars-svg").append("g")
        .attr("transform", `translate(${margin.xoffset},0)`)
      //
      var y0 = d3.scaleBand()
        .domain(group_names)
        .range([margin.top + 2, height - margin.bottom - 2])
        .paddingInner(0.1)
        .paddingOuter(0.0)
      //
      // var x = d3.scaleLog()
      var x = d3.scaleLinear()
        // .domain([d3.min(data, d => d.gene), d3.max(data, d => d.gene)]).nice()
        .domain([d3.min(bins, d => d.bin.min), d3.max(bins, d => d.bin.max)])
        .rangeRound([margin.left + 5, width - margin.right - 5])
      //
      var xAxis = g => g
        .attr("transform", `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).tickSizeOuter(0).ticks(5, "~g"))
        .call(g => g.select(".domain").remove())
        .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))
      // x-axis label right
      svg.append("text")
        .attr("x", (width / 2) - 30)
        .attr("y", height - margin.bottom / 5)
        .attr("alignment-baseline", "middle")
        .attr("font-family", "sans-serif")
        .attr("font-size", `${g_font_size}px`)
        .attr("font-weight", "bold")
        .text("Log2CPM")
      //
      var yAxis = g => g
        .attr("transform", `translate(${margin.left},0)`)
        .call(d3.axisLeft(y0).ticks(null, "s"))
        .call(g => g.select(".domain").remove())
        .call(g => g.selectAll(".tick text").attr("font-size", `${g_font_size}px`))
      // background stripes
      const bg = svg.append("g")
        .selectAll("g")
        .data(groups)
        .join("g")
      bg.append("rect")
        .attr("fill", (d, i) => i % 2 == 0 ? "#ffffff00" : "#eee")
        .attr("x", margin.left)
        .attr("width", width - margin.left - margin.right)
        .attr("y", d => y0(d[0]))
        .attr("height", y0.bandwidth())
      // 
      const g = svg.append("g")
        .selectAll("g")
        .data(bins)
        .join("g")
      // range lines
      g.append("path")
          .attr("transform",
            d => `translate(0,${y0(d[groupKey])})`
          )
          .attr("stroke", "currentColor")
          .attr("stroke-width", "0.5px")
          .attr("d", d => `
            M${x(d.bin.r0)},${y0.bandwidth() / 2}
            H${x(d.bin.r1)}
          `);
      g.append("g")
          .attr("transform",
            d => `translate(0,${y0(d[groupKey]) + y0.bandwidth() / 2})`
          )
        .selectAll("circle")
        .data(d => d.bin.outliers)
        .join("circle")
          .attr("cx", d => x(d))
          .attr("r", 1)
      // boxes
      let color_range = bins.length <= 8 ? pals.okabe : pals.mpn65
      let color = d3.scaleOrdinal()
        .domain(bins.map(d => d[groupKey]))
        .range(color_range)
      //
      g.append("path")
          .attr("fill", "#333")
          .attr("transform",
            d => `translate(0,${y0(d[groupKey])})`
          )
          // .attr("fill", d => pals.okabe[0])
          .attr("fill", d => color(d[groupKey]))
          .attr("stroke", "#000000")
          .attr("stroke-width", "0.5px")
          .attr("d", d => `
            M${x(d.bin.q3)},0
            V${y0.bandwidth()}
            H${x(d.bin.q1)}
            V${0}
            Z
          `);
      // median lines
      g.append("path")
          .attr("stroke", "currentColor")
          .attr("stroke-width", 2)
          .attr("transform",
            d => `translate(0,${y0(d[groupKey])})`
          )
          .attr("d", d => `
            M${x(d.bin.q2)},0
            V${y0.bandwidth()}
          `);
      //
      let panelBorder = svg.append("rect")
        .attr("x", margin.left)
        .attr("y", margin.top)
        .attr("height", height - margin.top - margin.bottom)
        .attr("width", width - margin.left - margin.right)
        .style("stroke", "black")
        .style("fill", "none")
        .style("stroke-width", "0.5px");
      //
      svg.append("text")
        .attr("x", margin.left)
        .attr("y", margin.top - 5)
        .attr("font-family", "sans-serif")
        .attr("font-size", `${g_font_size}px`)
        // .attr("font-style", "italic")
        .html(`Mean <tspan font-style="italic">${geneSym}</tspan> by ${aggKey}`)
      // draw
      svg.append("g")
          .call(xAxis);
      // svg.append("g")
      //     .call(yAxis);
    }

    // Use d3 to draw a binhex plot
    return {
      "drawGene": drawGene4,
      "drawMetaHex": drawMetaHex,
      "drawMetaBoxplot": draw_meta_boxplot,
      "drawMetaHeatmap": draw_meta_heatmap,
      "drawBars": drawBars,
      "drawGeneBars": draw_gene_bars,
      "drawGeneBoxplot": draw_gene_boxplot
    }
  }();

  function deparam(params){
  /* https://github.com/jupiterjs/jquerymx/blob/master/lang/string/deparam/deparam.js */
    if(! params || ! paramTest.test(params) ) {
      return {};
    }

    var data = {},
      pairs = params.split('&'),
      current;

    for(var i=0; i < pairs.length; i++){
      current = data;
      var pair = pairs[i].split('=');

      // if we find foo=1+1=2
      if(pair.length !== 2) {
        pair = [pair[0], pair.slice(1).join("=")];
      }

      var key = decodeURIComponent(pair[0].replace(plus, " ")),
        value = decodeURIComponent(pair[1].replace(plus, " ")),
        parts = key.match(keyBreaker);

      for ( var j = 0; j < parts.length - 1; j++ ) {
        var part = parts[j];
        if (!current[part] ) {
          // if what we are pointing to looks like an array
          current[part] = digitTest.test(parts[j+1]) || parts[j+1] === "[]" ? [] : {};
        }
        current = current[part];
      }
      var lastPart = parts[parts.length - 1];
      if(lastPart === "[]"){
        current.push(value);
      }else{
        current[lastPart] = value;
      }
    }
    return data;
  }

  function getVar(name, defVal) {
    /* get query variable from current URL or default value if undefined */
    var myUrl = window.location.href;
    const val = new URL(location.href).searchParams.get(name)
    if (val) {
      return val
    }
    return defVal
    // myUrl = myUrl.replace("#", "")
    // var urlParts = myUrl.split("?")
    // var queryStr = urlParts[1]
    // var varDict = deparam(queryStr) // parse key=val&... string to object
    // if (varDict[name] === undefined) {
    //   return defVal
    // }
    // return varDict[name]
  }

  function getDatasetNameFromUrl() {
    /* search for the "ds" parameter or a DNS hostname that indicates the dataset */
    // if ds=xxx was found in the URL, load the respective dataset
    var datasetName = getVar("ds");
    if (datasetName === undefined) {
        datasetName = ""
    }
    return datasetName;
  }

  function main(rootMd5) {
    console.log("mybrowser()");
    // var datasetName = 'Smillie2019';
    var datasetName = getDatasetNameFromUrl()
    state.ds = datasetName
    db = new CbDbFile(datasetName);
    db.loadConfig(onConfigLoaded, rootMd5);
  }

  function colorByDefaultField(onDone) {
    /* get the default coloring field from the config or the URL and start coloring by it.
     * Call onDone() when done. */
    var colorByInfo = db.getDefaultColorField();
    var colorType = colorByInfo[0];
    var colorBy = colorByInfo[1];

    if (g_loadGeneAndColor == 0) {
      renderer.drawMetaHex(colorBy)
    }

    if (g_donorField in g_mydata[0] && g_healthField in g_mydata[0]) {
      // renderer.drawMetaBoxplot(colorBy, state.gene_groupby)
      renderer.drawMetaBoxplot(colorBy, g_healthField)
      //renderer.drawMetaHeatmap(colorBy, g_healthField)
    }

    if (getVar("gene") !== undefined) {
      colorType = "gene";
      colorBy = getVar("gene");
      state.gene_groupby = getVar("groupby") || "none"
      loadGeneAndColor(colorBy) //"HMGB2")
    }

    // // allow to override coloring by URL args
    // else if (getVar("meta")!==undefined) {
    //   colorType = "meta";
    //   colorBy = getVar("meta");
    // }
    // if (colorType === "meta") {
    //   // update the meta field combo box
    //   var fieldIdx  = db.fieldNameToIndex(onlyAlphaNum(colorBy));
    //   if (fieldIdx===null) {
    //     alert("Default coloring is configured to be on field "+fieldName+
    //      " but cannot find a field with this name, using field 1 instead.");
    //     fieldIdx = 1;
    //   }
    //   $('#tpMetaCombo').val(fieldIdx).trigger('chosen:updated');
    //   $('#tpMetaBox_'+db.getMetaFields()[fieldIdx].name).addClass('tpMetaSelect');
    // }
    // else {
    //   loadGeneAndColor(colorBy, onDone);
    // }
  }

  var digitTest = /^\d+$/,
    keyBreaker = /([^\[\]]+)|(\[\])/g,
    plus = /\+/g,
    paramTest = /([^?#]*)(#.*)?$/;

  function deparam(params){
  /* https://github.com/jupiterjs/jquerymx/blob/master/lang/string/deparam/deparam.js */
    if(! params || ! paramTest.test(params) ) {
      return {};
    }

    var data = {},
      pairs = params.split('&'),
      current;

    for(var i=0; i < pairs.length; i++){
      current = data;
      var pair = pairs[i].split('=');

      // if we find foo=1+1=2
      if(pair.length !== 2) {
        pair = [pair[0], pair.slice(1).join("=")];
      }

      var key = decodeURIComponent(pair[0].replace(plus, " ")),
        value = decodeURIComponent(pair[1].replace(plus, " ")),
        parts = key.match(keyBreaker);

      for ( var j = 0; j < parts.length - 1; j++ ) {
        var part = parts[j];
        if (!current[part] ) {
          // if what we are pointing to looks like an array
          current[part] = digitTest.test(parts[j+1]) || parts[j+1] === "[]" ? [] : {};
        }
        current = current[part];
      }
      var lastPart = parts[parts.length - 1];
      if(lastPart === "[]"){
        current.push(value);
      }else{
        current[lastPart] = value;
      }
    }
    return data;
  }

  function changeUrl(vars, oldVars) {
  /* push the variables (object) into the history as the current URL. key=null deletes a variable. */
     // first get the current variables from the URL of the window
     var myUrl = window.location.href;
     myUrl = myUrl.replace("#", "");
     var urlParts = myUrl.split("?");
     var baseUrl = urlParts[0];

     let urlVars;
     if (oldVars===undefined) {
       var queryStr = urlParts[1];
       urlVars = deparam(queryStr); // parse key=val&... string to object
     } else {
       urlVars = oldVars;
     }

     // overwrite everthing that we got
     for (var key in vars) {
       var val = vars[key];
       if (val===null || val in urlVars)
         delete urlVars[key];
       else
         urlVars[key] = val;
     }

     var argStr = jQuery.param(urlVars); // convert to query-like string

     var dsName = "noname"
     if (db !== null) {
       dsName = db.getName()
     }

     if (argStr.length > 1000) {
       warn("Cannot save current changes to the URL, the URL would be too long. "+
         "You can try to shorten some cluster labels to work around the problem temporarily. "+
         "But please contact us at cells@ucsc.edu and tell us about the error. Thanks!");
     } else {
       if (vars.ds) {
         // window.location.href = baseUrl+"?"+argStr
         history.replaceState({}, dsName, baseUrl+"?"+argStr);
         state.ds = vars.ds
         state.gene_loaded = false
				 db = new CbDbFile(vars.ds)
				 db.loadConfig(onConfigLoaded, dsName)
       } else {
         history.replaceState({}, dsName, baseUrl+"?"+argStr);
       }
     }
  }

  function gotFirstCoords(coords, info, clusterMids) {
    /* XX very ugly way to implement promises. Need a better approach one day. */
    // gotCoords(coords, info, clusterMids);
    var n_cells = coords.length / 2;
    // console.log(n_cells + ' cells');
    // console.log(coords);
    // console.log(info);
    // console.log(clusterMids);

    // if (g_mydata === null) {
    //   for (var i = 0; i < coords.length; i += 2) {
    //     g_mydata.push({x: coords[i], y: coords[i + 1]});
    //   }
    // } else {
    //   for (var i = 0; i < n_cells; i++) {
    //     g_mydata[i].x = coords[i * 2];
    //     g_mydata[i].y = coords[i * 2 + 1];
    //   }
    // }

    g_mydata = []
    for (var i = 0; i < coords.length; i += 2) {
      g_mydata.push({ x: coords[i], y: coords[i + 1] });
    }

    function gotMetaArr(metaArr, metaInfo, funcVal) {
      var fieldName = metaInfo.name;
      for (var i = 0; i < metaArr.length; i++) {
        g_mydata[i][fieldName] = metaInfo.valCounts[metaArr[i]][0];
      }
      colorByDefaultField()
    }

    // Set the title of the page
    setTitle(db.conf.shortLabel)

		setupDataChooser()

    var metaInfo = db.findMetaInfo(db.conf.clusterField)
    db.loadMetaVec(metaInfo, gotMetaArr, onProgressConsole)
    
    var meta_fields = db.conf.metaFields
      .filter(
        d => d.type == "enum" &&
        (d.valCounts.length <= 8 && d.name != "custom") |
        (d.name == g_donorField)
      )
    for (const meta_field of meta_fields) {
      db.loadMetaVec(db.findMetaInfo(meta_field.name), gotMetaArr, onProgressConsole)
    }

    console.log("updating gene plots")
    console.log({state})
    console.log({g_mydata})
    renderer.drawGene(g_mydata, state.gene, state.gene_groupby)
    renderer.drawGeneBars(state.gene, state.gene_groupby)
    renderer.drawGeneBoxplot(state.gene, state.gene_groupby)
    renderer.drawMetaHex(state.meta_colorby)
    renderer.drawMetaBoxplot(state.meta_colorby, g_healthField)
    // getGeneInfo(state.gene)


  } // gotFirstCoords()

  function setTitle(text) {
    // $("#mytitle").text(text)
  }

	function setupDataChooser() {
    var tissue_dirs = [
      {'id': 'a12_4_4_t4_cd8_1_2', 'label': 'CD8 T cells'},
      {'id': 'a12_4_4_t4_cd4_2_2', 'label': 'CD4 T cells'},
      {'id': 'a12_4_4_m3_2', 'label': 'Myeloid cells'},
      {'id': 'a12_4_4_b5_1_3', 'label': 'B cells'},
      {'id': 'n3_2', 'label': 'Epithelial nuclei'},
    ]
    var blood_dirs = [
      {'id': 'blood2_tcell5_cd8_5', 'label': 'CD8 T cells'},
      {'id': 'blood2_tcell5_cd4_5', 'label': 'CD4 T cells'},
      {'id': 'blood2_myeloid5', 'label': 'Myeloid cells'},
      {'id': 'blood2_bcell5', 'label': 'B cells'}
    ]

    var dsName = db.getName()

    var htmls = [
		  '<div class="form-check form-check-inline"><b>Tissue</b></div>'
    ]
    for (let dir of tissue_dirs) {
      var x = `<div class="form-check form-check-inline">
      <input class="form-check-input" type="radio" name="input-choose-tissue" id="${dir.id}" value="${dir.id}" ${dsName == dir.id ? 'checked' : ''}>
      <label class="form-check-label ${dsName == dir.id ? 'font-weight-bold' : ''}" for="${dir.id}">${dir.label}</label>
      </div>`
      htmls.push(x)
    }
    $("#choose-tissue").html(htmls.join(""))

    var htmls = [
		  '<div class="form-check form-check-inline"><b>Blood</b></div>'
    ]
    for (let dir of blood_dirs) {
      var x = `<div class="form-check form-check-inline">
      <input class="form-check-input" type="radio" name="input-choose-tissue" id="${dir.id}" value="${dir.id}" ${dsName == dir.id ? 'checked' : ''}>
      <label class="form-check-label ${dsName == dir.id ? 'font-weight-bold' : ''}" for="${dir.id}">${dir.label}</label>
      </div>`
      htmls.push(x)
    }
    $("#choose-blood").html(htmls.join(""))

    $("[name='input-choose-tissue']").change(
      ({target}) => {
        let ds = target.getAttribute('value')
        // changeUrl({"gene":state.gene, "ds":ds, "groupby":state.gene_groupby})
        changeUrl({"ds": ds})
      }
    )
	}

  function makeLabelRenames(metaInfo) {
    /* return an obj with old cluster name -> new cluster name */
    var valCounts = metaInfo.valCounts;
    if (valCounts===undefined) { // 'int' and 'float' types do not have their values counted yet
      // this doesn't work because the values are not loaded yet, requires moving this call to 
      // later
      //metaInfo.valCounts = countValues(metaInfo.arr);
      alert("cannot label on numeric fields, please use the enumFields option in cellbrowser.conf");
    }
    var newLabels = metaInfo.ui.shortLabels;

    var oldToNew = {};
    for (var i = 0; i < valCounts.length; i++) {
      var oldLabel = valCounts[i][0];
      var newLabel = newLabels[i];
      oldToNew[oldLabel] = newLabel;
    }
    return oldToNew;
  }

  function onConfigLoaded(datasetName) { 
    // this is a collection if it does not have any field information
    if (!db.conf.metaFields) {
      console.log("could not find db.conf.metaFields")
      return
    }

    let binData = localStorage.getItem(db.name+"|custom");
    if (binData) {
      let jsonStr = LZString.decompress(binData);
      let customMeta = JSON.parse(jsonStr);
      db.conf.metaFields.unshift(customMeta);
    }

    db.loadCoords(0, gotFirstCoords, onProgressConsole);
    // Show the table of gene markers for each cluster
    // onClusterNameClick("TA 1")
    var metaInfo = db.findMetaInfo(db.conf.clusterField)
    console.log("onConfigLoaded()")
    console.log({metaInfo})
    var clusterNames = metaInfo.valCounts.map(d => d[0]).slice().sort()
    buildMetaControls()
    buildMarkerTables(() => {
      console.log(`finished buildMarkerTables(), ${clusterNames[0]}`)
      onClusterNameClick(0, 0, clusterNames[0])
    })
  }

  function onGeneChange(ev) {
    /* user changed the gene in the combobox */
    var geneSym = ev.target.value;
    if (geneSym === "") {
      return; // do nothing if user just deleted the current gene
    }
    loadGeneAndColor(geneSym);
  }

  function geneComboSearch(query, callback) {
    /* called when the user types something into the gene box, returns matching gene symbols */
    if (!query.length) {
      return callback()
    }
    console.log("LOG geneComboSearch")
    db.searchGenes(query.toLowerCase(), callback)
  }

  function onProgressConsole(ev) {
    console.log(ev);
  }

  function loadGeneAndColor(geneSym, onDone) {
    if (state.gene == geneSym && state.gene_loaded) {
      return
    }
    state.gene = geneSym
    // document.getElementById("link-gene-contrasts").href = `/ircolitis/gene-contrasts/?gene=${state.gene}`
    document.getElementById("link-gene-contrasts").href = `/ircolitis/gene-contrasts/`
    recent_genes.unshift(geneSym)
    recent_genes = recent_genes.filter((x, i, a) => a.indexOf(x) === i).slice(0, 10)
    $("#recent-genes").html(
      "Recent: " + recent_genes.map(x =>
        `<i><a data-gene="${x}" class="text-primary recent-gene-link" style="cursor: pointer;">${x}</a></i>`
      ).join(", ")
    )
    $(".recent-gene-link").unbind("click")
    $(".recent-gene-link").on("click", onMarkerGeneClick);

    ++g_loadGeneAndColor
    console.log('loadGeneAndColor() call ' + g_loadGeneAndColor)
    /* color by a gene, load the array into the renderer and call onDone or just redraw */
    if (onDone===undefined) {
      onDone = function() {
        // if (state.gene != geneSym) {
        //   state.gene = geneSym
          state.gene_loaded = true
          renderer.drawGene(g_mydata, geneSym, state.gene_groupby)
          renderer.drawGeneBars(geneSym, state.gene_groupby)
          renderer.drawGeneBoxplot(geneSym, state.gene_groupby)
					getGeneInfo(geneSym)
        // }
      }
    }

    function gotGeneVec(exprArr, decArr, geneSym, geneDesc, binInfo) {

      g_exprArr = exprArr;
      g_decArr = decArr;
      g_geneSym = geneSym;
      g_geneDesc = geneDesc;
      g_binInfo = binInfo;

      if (g_mydata === null) {
        for (var i = 0; i < exprArr.length; i++) {
          g_mydata.push({gene: exprArr[i]});
        }
      } else {
        for (var i = 0; i < exprArr.length; i++) {
          g_mydata[i]['gene'] = exprArr[i];
        }
      }
      
      /* called when the expression vector has been loaded and binning is done */
      if (decArr===null)
        return;
      console.log("Received expression vector, gene "+geneSym+", geneId "+geneDesc)
      // _dump(binInfo);
      // makeLegendExpr(geneSym, geneDesc, binInfo, exprArr, decArr);

      // renderer.setColors(legendGetColors(gLegend.rows));
      // renderer.setColorArr(decArr);
      // buildLegendBar();
      onDone()

      // update the gene combo box
      // selectizeSetValue("#tpGeneCombo", geneSym);

      // // update the "recent genes" div
      // for (var i = 0; i < gRecentGenes.length; i++) {
      //   // remove previous gene entry with the same symbol
      //   if (gRecentGenes[i][0]===geneSym) {
      //     gRecentGenes.splice(i, 1);
      //     break;
      //   }
      // }
      // gRecentGenes.unshift([geneSym, geneDesc]); // insert at position 0
      // gRecentGenes = gRecentGenes.slice(0, 9); // keep only nine last
    }

    changeUrl({"gene":geneSym, "groupby":state.gene_groupby});
    console.log("Loading gene expression vector for "+geneSym);

    db.loadExprAndDiscretize(geneSym, gotGeneVec, onProgressConsole);
  }

  function buildMetaControls() {
    var htmls = []

    // Checkbox for showing or hiding cluster labels on the UMAP
    htmls.push(`<div class="col-12 col-md-5">
      <input type="checkbox" id="show-cluster-labels" checked>
      <label for="show-cluster-labels">Cluster labels</label>
    </div>`)

    // Dropdown select for what metadata variable to use as color on the UMAP
    var div_colorby = [
      `<div class="col-12 col-md-5">
        <div class="form-group">
          <label for="meta-colorby">Color by</label>
          <select class="form-control" id="meta-colorby">
      `
    ]
    for (var i = 0; i < db.conf.metaFields.length; i++) {
      var fieldLabel = db.conf.metaFields[i].label
      var fieldName = db.conf.metaFields[i].name
      if (fieldName.startsWith("custom") || fieldName.toLowerCase().startsWith("cell")) {
        continue
      }
      if (fieldName === db.conf.clusterField) {
        div_colorby.push(`<option selected="selected" value="${fieldName}">${fieldLabel}</option>`)
      } else {
        div_colorby.push(`<option value="${fieldName}">${fieldLabel}</option>`)
      }
    }
    div_colorby.push(`
          </select>
        </div>
      </div>
    `)
    div_colorby = div_colorby.join("")
    htmls.push(div_colorby)

    $("#meta-controls").html(htmls.join(""))

    htmls = []
    var div_groupby = [
      `<div class="col-12 col-lg-4">
        <div class="form-group">
          <label for="gene-groupby">Group by</label>
          <select class="form-control" id="gene-groupby" onchange="mybrowser.onGeneGroupbyClick()">
      `
    ]
    var meta_fields = db.conf.metaFields
      .filter(
          d => d.type == "enum" && d.valCounts.length <= 8 &&
               d.name != "custom" && !d.name.startsWith("leiden")
      )
    // div_groupby.push(`<option value="none">None</option>`)
    for (let meta_field of [...meta_fields, {name: "none", label: "None"}]) {
      div_groupby.push(`<option value="${meta_field.name}" ${meta_field.name == getVar("groupby") ? "selected" : ""}>${meta_field.label}</option>`)
    }
    div_groupby.push(`
          </select>
        </div>
      </div>
    `)
    div_groupby = div_groupby.join("")

    // Gene expression controls
    htmls.push(`
      <div class="col-12 col-lg-8">
        <label for="gene-search">Gene</label>
        <select style="max-width:400px" id="gene-search" placeholder="search for a gene..." class="tpCombo"></select>
        <div class="" id="recent-genes"></div>
      </div>
      ${div_groupby}
    `); // Genes

    $("#gene-controls").html(htmls.join(""))

    let gene_search = $('#gene-search').selectize({
      maxItems: 1,
      valueField : 'id',
      labelField : 'text',
      searchField : 'text',
      closeAfterSelect: false,
      load : geneComboSearch
    });
    gene_search.on("change", onGeneChange);
    $("#show-cluster-labels").unbind("change")
    $("#show-cluster-labels").change(function() {$(".cluster_label").toggle()})
    $("#meta-colorby").unbind("change")
    $("#meta-colorby").on("change", on_meta_colorby)
  }

  function buildMarkerTables(onDone = () => null) {

    // e.g. "1", "2", "1-vs-2", ...
    var tab_names = db.conf.markers.map(d => sanitizeName(d.shortLabel.split("|")[0]))

    // e.g. "ova", "ava", "con"
    var tab_types = db.conf.markers.map(d => d.shortLabel.split("|")[1])

    // Choose the resolution with the matching cluster names
    var this_tab = db.conf.clusterField
    if (db.conf.markers[0].shortLabel.includes("|")) {
      if (this_tab === "cluster") {
        // this_tab = tab_names.slice(-1)[0]
        var metaInfo = db.findMetaInfo(this_tab)
        var clusterNames = metaInfo.valCounts.map(d => d[0])
        var a = new Set(clusterNames)
        for (let i = 0; i < tab_names.length; i++) {
          var b = new Set(db.conf.markers[i].shortLabel.split("|")[3].split("\t"))
          let intersection = [...a].filter(x => b.has(x))
          if (intersection.length == a.size) {
            this_tab = tab_names[i]
            break
          }
        }
      }
    }

    let tab_ova = -1
    let tab_ava = -1
    for (let i = 0; i < tab_names.length; i++) {
      if (this_tab == tab_names[i]) {
        if (tab_types[i] == "ova") {
          tab_ova = i
        } else if (tab_types[i] == "ava") {
          tab_ava = i
        }
      }
    }
    console.log(`buildMarkerTables() for ${this_tab} tab_ova=${tab_ova} tab_ava=${tab_ava}`)

    var htmls = []

    if (tab_ova == -1) {

      var metaInfo = db.findMetaInfo(db.conf.clusterField)
      var clusterNames = metaInfo.valCounts.map(d => d[0]).slice().sort((a, b) => a.localeCompare(b, navigator.languages[0] || navigator.language, {numeric: true, ignorePunctuation: true}))

      // htmls.push(`<div><h4>Clusters</h4></div>`)

      htmls.push(htmlTabs(0, clusterNames))

      // htmls.push(`<nav> <div class="nav nav-pills" id="nav-cluster-tables" role="tablist"> `)
      // for (var clusterIndex = 0; clusterIndex < clusterNames.length; clusterIndex++) {
      //   var clusterName = clusterNames[clusterIndex]
      //   htmls.push(`<a class="nav-item nav-link ${clusterIndex == 0 ? 'active' : ''}" id="nav-cluster-table-${clusterIndex}" data-toggle="tab" href="#cluster-table-${clusterIndex}" role="tab" aria-controls="cluster-table-${clusterIndex}" aria-selected="${clusterIndex == 0 ? true : false}" onclick="mybrowser.onClusterNameClick(${clusterIndex}, '${clusterName}')">${clusterName}</a>`)
      // }
      // htmls.push(`</div> </nav>`)
      // htmls.push("<div id='tpPaneHeader' style='padding:0.4em 1em'>");
      // htmls.push("</div>");
      // htmls.push("<div class='tab-content' id='cluster-tables'>");
      // for (var clusterIndex = 0; clusterIndex < clusterNames.length; clusterIndex++) {
      //   htmls.push(`<div class="tab-pane ${clusterIndex == 0 ? 'show active' : ''}" id="cluster-table-${clusterIndex}" role="tabpanel" aria-labelledby="nav-cluster-table-${clusterIndex}">`)
      //   htmls.push("</div>"); // tab-pane
      // }
      // htmls.push("</div>"); // tab-content
      $("#mytable").html(htmls.join(""));

    } else if (tab_ova != -1) {

      // htmls.push(`<div><h4>Genes</h4></div>`)

      var ova_names = db.conf.markers[tab_ova].shortLabel.split("|")[3].split("\t")
      var ova_tabs = htmlTabs(tab_ova, ova_names)

      var ava_tabs = ""
      if (tab_ava != -1) {
        var ava_names = db.conf.markers[tab_ava].shortLabel.split("|")[3].split("\t")
        ava_tabs = htmlTabs(tab_ava, ava_names)
      }

      htmls.push(`
      <ul class="nav nav-tabs" id="myTab" role="tablist">
        <li class="nav-item" role="presentation">
          <a class="nav-link active" id="ova-tab" data-toggle="tab" href="#ova" role="tab" aria-controls="ova" aria-selected="true">One vs All</a>
        </li>
      `)
      if (tab_ava != -1) {
        htmls.push(`
          <li class="nav-item" role="presentation">
            <a class="nav-link" id="ava-tab" data-toggle="tab" href="#ava" role="tab" aria-controls="ava" aria-selected="false">All vs All</a>
          </li>`)
      }
      htmls.push(`
      </ul>
      <div class="tab-content">
        <div class="tab-pane active" id="ova" role="tabpanel" aria-labelledby="ova-tab">${ova_tabs}</div>
      `)
      if (tab_ava != -1) {
        htmls.push(`
          <div class="tab-pane" id="ava" role="tabpanel" aria-labelledby="ava-tab">${ava_tabs}</div>
        `)
      }
      htmls.push(`</div>`)
      
      $("#mytable").html(htmls.join(""));
    
    }

    onDone()
  } // buildMarkerTables

  function htmlTabs(tabId, clusterNames) {
    var htmls = []
    var style_overflow = ""
    if (clusterNames.length > 30) {
      style_overflow = "style='max-height:6rem;overflow-y:scroll;'"
    }
    htmls.push(`<nav> <div ${style_overflow} class="nav nav-pills" id="nav-cluster-tables-${tabId}" role="tablist"> `)
    for (var clusterIndex = 0; clusterIndex < clusterNames.length; clusterIndex++) {
      var clusterName = clusterNames[clusterIndex]
      htmls.push(`<a class="cluster-tab nav-item nav-link ${clusterIndex == 0 ? 'active' : ''}" id="nav-cluster-table-${tabId}-${clusterIndex}" data-toggle="tab" href="#cluster-table-${tabId}-${clusterIndex}" role="tab" aria-controls="cluster-table-${tabId}-${clusterIndex}" aria-selected="${clusterIndex == 0 ? true : false}" onclick="mybrowser.onClusterNameClick(${tabId}, ${clusterIndex}, '${clusterName}')">${clusterName}</a>`)
    }
    htmls.push(`</div> </nav>`)

    htmls.push("<div id='tpPaneHeader' style='padding:0.4em 1em'>");
    htmls.push("</div>");

    htmls.push("<div class='tab-content' id='cluster-tables'>");
    for (var clusterIndex = 0; clusterIndex < clusterNames.length; clusterIndex++) {
      htmls.push(`<div class="tab-pane ${clusterIndex == 0 ? 'show active' : ''}" id="cluster-table-${tabId}-${clusterIndex}" role="tabpanel" aria-labelledby="nav-cluster-table-${tabId}-${clusterIndex}">`)
      htmls.push("</div>"); // tab-pane
    }
    htmls.push("</div>"); // tab-content
    return htmls.join("")
  }

  function onGeneGroupbyClick() {
    var x = document.getElementById("gene-groupby").value
    // if (x != "none") {
    //   $("#gene-legend").html("").css("min-height", "120px")
    // } else {
    //   $("#gene-legend").html("").css("min-height", "0px")
    // }
    state.gene_groupby = x
    changeUrl({"groupby":state.gene_groupby})
    renderer.drawGene(g_mydata, state.gene, state.gene_groupby)
    renderer.drawGeneBars(state.gene, state.gene_groupby)
    renderer.drawGeneBoxplot(state.gene, state.gene_groupby)
    console.log("changed to " + x)
  }

  function onClusterNameClick(tabId, clusterIndex, clusterName) {

    console.log("onClusterNameClick() - building marker genes window for " + clusterName)
    var htmls = []
    var divName = `tabs-${tabId}-cluster-${clusterIndex}`
    var sanName = sanitizeName(clusterName.replaceAll('-', 'Minus'))
    var markerTsvUrl = cbUtil.joinPaths([
      db.name, "markers", db.conf.markers[tabId].name, `${sanName}.tsv.gz`
    ])
    $(".my-gene-table").remove()
    htmls.push(`<div class="my-gene-table" id="${divName}" style="height:500px;">Loading...</div>`)
    $(`#cluster-table-${tabId}-${clusterIndex}`).html(htmls.join(""))
    loadClusterTsv(markerTsvUrl, loadMarkersFromTsv2, divName, clusterName, clusterIndex)

    // var metaInfo = db.findMetaInfo(db.conf.clusterField)

    // var tabInfo = db.conf.markers; // list with (label, subdirectory)
    // var htmls = []

    // // e.g. "1", "2", "1-vs-2", ...
    // var tab_names = db.conf.markers.map(d => sanitizeName(d.shortLabel.split("|")[0]))

    // // e.g. "ova", "ava", "con"
    // var tab_types = db.conf.markers.map(d => d.shortLabel.split("|")[1])

    // // Choose the last name (highest resolution) if the name is "cluster"
    // var this_tab = db.conf.clusterField
    // if (this_tab === "cluster") {
    //   this_tab = tab_names.slice(-1)[0]
    // }

    // let tab_ova = -1
    // let tab_ava = -1

    // // Find which tables correspond to this resolution
    // for (let i = 0; i < tabInfo.length; i++) {
    //   var sanTab = sanitizeName(tabInfo[i].shortLabel.split("|")[0])
    //   if (this_tab == tab_names[i]) {
    //     if (tab_types[i] == "ova") {
    //       tab_ova = i
    //     } else if (tab_types[i] == "ava") {
    //       tab_ava = i
    //     }
    //   }
    // }

    // if (tab_ova != -1) {
    //   var divName = `tabs-${tabId}-cluster-${clusterIndex}`
    //   var sanName = sanitizeName(clusterName)
    //   var markerTsvUrl = cbUtil.joinPaths([
    //     db.name, "markers", tabInfo[tab_ova].name, `${sanName}.tsv.gz`
    //   ])
    //   $(".my-gene-table").remove()
    //   htmls.push(`<div class="my-gene-table" id="${divName}" style="height:500px;">Loading...</div>`)
    //   $(`#cluster-table-${tabId}-${clusterIndex}`).html(htmls.join(""))
    //   loadClusterTsv(markerTsvUrl, loadMarkersFromTsv2, divName, clusterName, clusterIndex)
    // }

  }

  // ---------

  var DEBUG = true;

  function _dump(o) {
  /* for debugging */
    console.log(JSON.stringify(o));
  }

  function formatString (str) {
    /* Stackoverflow code https://stackoverflow.com/a/18234317/233871 */
    /* "a{0}bcd{1}ef".formatUnicorn("foo", "bar"); // yields "aFOObcdBARef" */
      if (arguments.length) {
        var t = typeof arguments[0];
        var key;
        var args = ("string" === t || "number" === t) ?
          Array.prototype.slice.call(arguments)
          : arguments[0];
  
        for (key in args) {
          str = str.replace(new RegExp("\\{" + key + "\\}", "gi"), args[key]);
        }
      }
      return str;
    }

  function debug(msg, args) {
    if (DEBUG) {
      console.log(formatString(msg, args));
    }
  }

  function warn(msg) {
    alert(msg);
  }

  function cloneObj(d) {
  /* returns a copy of an object, wasteful */
    // see http://stackoverflow.com/questions/122102/what-is-the-most-efficient-way-to-deep-clone-an-object-in-javascript
    return JSON.parse(JSON.stringify(d));
  }

  function cloneArray(a) {
  /* returns a copy of an array */
    return a.slice();
  }

  function copyNonNull(srcArr, trgArr) {
  /* copy non-null values to trgArr */
    if (srcArr.length!==trgArr.length)
      alert("warning - copyNonNull - target and source array have different sizes.");

    for (var i = 0; i < srcArr.length; i++) {
      if (srcArr[i]!==null)
        trgArr[i] = srcArr[i];
    }
    return trgArr;
  }

  function isEmpty(obj) {
    for(var key in obj) {
      if(obj.hasOwnProperty(key))
        return false;
    }
    return true;
  }

  function allEmpty(arr) {
    /* return true if all members of array are white space only strings */
    var newArr = arr.filter(function(str) { return /\S/.test(str); });
    return (newArr.length===0);
  }

  function copyNonEmpty(srcArr, trgArr) {
  /* copy from src to target array if value is not "". Just return trgArr is srcArr is null or lengths don't match.  */
    if (!srcArr || (srcArr.length!==trgArr.length))
      return trgArr;

    for (var i = 0; i < srcArr.length; i++) {
      if (srcArr[i]!=="")
        trgArr[i] = srcArr[i];
    }
    return trgArr;
  }

  function keys(o) {
  /* return all keys of object as an array */
    var allKeys = [];
    for(var k in o) allKeys.push(k);
    return allKeys;
  }

  function capitalize(s) {
    return s[0].toUpperCase() + s.slice(1);
  }

  
  function findMetaValIndex(metaInfo, value) {
    /* return the index of the value of an enum meta field */
    var valCounts = metaInfo.valCounts;
    for (var valIdx = 0; valIdx < valCounts.length; valIdx++) {
      if (valCounts[valIdx][0]===value)
        return valIdx;
    }
  }

  function intersectArrays(arrList) {
    /* return the intersection of all arrays as an array. Non-IE11? */
    var smallSet = new Set(arrList[0]);
    for (var i=1; i < arrList.length; i++) {
      var otherSet = new Set(arrList[i]);
      smallSet = new Set([...smallSet].filter(x => otherSet.has(x)));
    }
    var newArr = Array.from(smallSet);
    // alternative without spread:
    //function intersection(setA, setB) {
      //  var _intersection = new Set();
      //  for (var elem of setB) {
      //    if (setA.has(elem)) {
      //      _intersection.add(elem);
      //    }
      //  }
      //  return _intersection;
    //}
    return newArr;
  }

  function prettyNumber(/*int*/ count) /*str*/ {
    /* convert a number to a shorter string, e.g. 1200 -> 1.2k, 1200000 -> 1.2M, etc */
    var f = count;

    if (count>1000000) {
      f = (count / 1000000);
      return f.toFixed(1)+"M";
    }
    if (count>10000) {
      f = (count / 1000);
      return f.toFixed(0)+"k";
    }
    if (count>1000) {
      f = (count / 1000);
      return f.toFixed(1)+"k";
    }

    return f;
  }

  function sanitizeName(name) {
    /* ported from cellbrowser.py: remove non-alpha, allow underscores */
    var newName = name.replace(/[^a-zA-Z_0-9]/g, "");
    return newName;
  }

  function onlyAlphaNum(name) {
    /* only allow alphanumeric characters */
    var newName = name.replace(/[^a-zA-Z0-9+]/g, "");
    return newName;
  }

  function loadClusterTsv(fullUrl, func, divName, clusterName, clusterIndex) {
    console.log(`loadClusterTsv() clusterName=${clusterName} markerTsvUrl=${fullUrl}`)
  /* load a tsv file relative to baseUrl and call a function when done */
    function conversionDone(data) {
      Papa.parse(data, {
          complete: function(results, localFile) {
                func(results, localFile, divName, clusterName, clusterIndex);
              },
          error: function(err, file) {
                if (divName!==undefined)
                  alert("could not load "+fullUrl);
              }
          });
    }

    function onTsvLoadDone(res) {
      var data = res.target.response;
      if (res.target.responseURL.endsWith(".gz")) {
        data = pako.ungzip(data);
        //data = String.fromCharCode.apply(null, data); // only good for short strings
        data = arrayBufferToString(data, conversionDone);
      }
      else
        conversionDone(data);
    }

    var req = new XMLHttpRequest();
    req.addEventListener("load", onTsvLoadDone);
    req.open('GET', fullUrl, true);
    req.responseType = "arraybuffer";
    req.send();
  }

  function loadMarkersFromTsv2(papaResults, url, divId, clusterName, clusterIndex) {
    /* construct a table from a marker tsv file and write as html to the DIV with divID */
    console.log(`got coordinate TSV rows, parsing... divId=${divId} clusterName=${clusterName} clusterIndex=${clusterIndex}`)
    var rows = papaResults.data;
    var headerRow = rows[0];
    var columns = []
    var format_gene = function(cell, params, onRendered) {
      // return `<i><a data-gene="${cell.getValue()}" class='text-primary tpLoadGeneLink' style="cursor: pointer;">${cell.getValue()}</a></i>`
      onRendered(function(){
        var el = $(cell.getElement())
        el.html(`<i><a data-gene="${cell.getValue()}" class='text-primary tpLoadGeneLink' style="cursor: pointer;">${cell.getValue()}</a></i>`)
        el.click(onMarkerGeneClick)
      })
    }
    var format_numeric = function(cell, params, onRendered) {
      var v = cell.getValue() 
      if (v) {
        if (v > 0 && v < 0.01) {
          return v.toExponential(1)
        }
        return v.toFixed(2)
      }
      return ''
    }
    const pval_labels = [
      'P_value', 'pval', 'pvalue', 'P.Value', 'adj.P.Val', 'p_val', 'pVal', 'Chisq_P', 'fdr', 'FDR'
    ]
    for (var col = 1; col < headerRow.length; col++) {
      var coldata = headerRow[col].split("|")
      var d = {
        title: coldata[0],
        field: coldata[0]
      }
      if (coldata[0] === "symbol") {
        d["formatter"] = format_gene
        d["headerFilter"] = true
      }
      if (coldata[1] === "float") {
        d["headerFilter"] = true
        d["headerFilterFunc"] = ">="
        if (coldata[0] == "pct_out" || pval_labels.indexOf(coldata[0]) != -1) {
          d["headerFilterFunc"] = "<="
        }
        d["formatter"] = format_numeric
      }
      columns.push(d)
    }
    var data = []
    for (var i = 1; i < rows.length; i++ ) {
      var d = {}
      for (var col = 1; col < headerRow.length; col++) {
        var coldata = headerRow[col].split("|")
        var val = rows[i][col]
        if (coldata[1] === "float") {
          // val = (+val).toFixed(3)
          val = +val
        }
        d[coldata[0]] = val
      }
      if (rows[i][0] !== "") {
        data.push(d)
      }
    }
    $(`#${divId}`).html("")
    if ($(`#${divId}`).length) {
      var table = new Tabulator(`#${divId}`, {
        selectable: false,
        // set height of table (in CSS or here), this enables the Virtual DOM and
        // improves render speed dramatically (can be any valid css height value)
        height: 500,
        data: data,
        // layout: "fitDataFill",
        layout: "fitColumns",
        // columnMaxWidth: 290,
        movableColumns: true,
        columns: columns,
        resizableColumns: false
      })
    }
    $(`#${divId}`).addClass("table-striped")
  }

  function onMarkerGeneClick(ev) {
    console.log("onMarkerGeneClick(ev)")
    console.log({ev})
    /* user clicks onto a gene in the table of the marker gene dialog window */
    var geneSym = ev.target.getAttribute("data-gene");
    console.log(`geneSym = ${geneSym}`)
    if (geneSym) {
      loadGeneAndColor(geneSym);
    }
  }

  //asdf
  function on_meta_colorby() {
    var x = document.getElementById("meta-colorby").value
    state.meta_colorby = x
    console.log("meta_colorby changed to " + x)
    var onDone = () => null
    if (
      db.conf.clusterField != x &&
      (x == "cluster" || x.startsWith("leiden"))
    ) {
      db.conf.clusterField = x
      onDone = function() {
        buildMarkerTables(() => {
          $('[aria-selected="true"]')[1].click()
        })
        renderer.drawGeneBars(state.gene, state.gene_groupby)
        renderer.drawGeneBoxplot(state.gene, state.gene_groupby)
        renderer.drawMetaBoxplot(state.meta_colorby, g_healthField)
      }
    }
    renderer.drawMetaHex(state.meta_colorby, onDone)
  }

  // function onMetaClick(ev) {
  //   /* user clicks onto a gene in the table of the marker gene dialog window */
  //   var fieldName = ev.target.getAttribute("data-name");
  //   renderer.drawMetaHex(fieldName);
  //   // renderer.drawMetaBoxplot(fieldName, state.gene_groupby)
  //   renderer.drawMetaBoxplot(fieldName, g_healthField)
  // }

  function buildGeneCombo(htmls, id, left, width) {
    /* Combobox that allows searching for genes */
    //htmls.push('<div class="tpToolBarItem" style="position:absolute;left:'+left+'px;top:'+toolBarComboTop+'px">');
    htmls.push('<div class="tpToolBarItem" style="padding-left: 3px">');
    htmls.push('<label style="display:block; margin-bottom:8px; padding-top: 8px;" for="'+id+'">Color by Gene</label>');
    htmls.push('<select style="width:'+width+'px" id="'+id+'" placeholder="search for gene..." class="tpCombo">');
    htmls.push('</select>');
    htmls.push('</div>');
    //htmls.push("<button>Multi-Gene</button>");
  }

  function arrayBufferToString(buf, callback) {
    /* https://stackoverflow.com/questions/8936984/uint8array-to-string-in-javascript */
     var bb = new Blob([new Uint8Array(buf)]);
     var f = new FileReader();
     f.onload = function(e) {
       callback(e.target.result);
     };
     f.readAsText(bb);
  }

  function linspace(start, stop, nsteps) {
    let delta = (stop - start) / (nsteps - 1)
    return d3.range(nsteps).map(function(i){return start + i * delta;})
  }

  /**
   * Uses canvas.measureText to compute and return the width of the given text
   * of given font in pixels.
   * 
   * @param {String} text The text to be rendered.
   * @param {String} font The css font descriptor that text is to be rendered
   * with (e.g. "bold 14px verdana").
   * 
   * @see https://stackoverflow.com/questions/118241/calculate-text-width-with-javascript/21015393#21015393
   */
  function getTextWidth(text, font) {
      // re-use canvas object for better performance
      var canvas = getTextWidth.canvas || (getTextWidth.canvas = document.createElement("canvas"));
      var context = canvas.getContext("2d");
      context.font = font;
      var metrics = context.measureText(text);
      return metrics.width;
  }

  function max_length(xs) {
    var retval = 0
    var ret_i = 0
    for (var i = 0; i < xs.length; i++) {
      retval = xs[i].length > retval ? xs[i].length : retval
    }
    return retval
  }

  function ramp(color, n = 256) {
    // const canvas = DOM.canvas(1, n);
    var canvas = document.createElement("canvas")
    canvas.width = 1
    canvas.height = n
    const context = canvas.getContext("2d");
    for (let i = 0; i < n; ++i) {
      context.fillStyle = color(i / (n - 1));
      context.fillRect(0, n-i, 1, 1);
    }
    return canvas;
  }

  function vertical_legend({
    svg,
    color,
    title,
    tickSize = 3,
    width = 26 + tickSize,
    height = 320,
    offsetLeft = 0,
    offsetTop = 0,
    marginTop = 10,
    marginRight = 10 + tickSize,
    marginBottom = 10,
    marginLeft = 5,
    fontSize = "12px",
    ticks = height / 64,
    tickFormat,
    tickValues
  } = {}) {

    // let tickAdjust = g => g.selectAll(".tick line")
      // .attr("x1", marginLeft - width + marginRight);
    let tickAdjust = g => g.selectAll(".tick text")
      .attr("font-size", fontSize)
      .attr("font-family", "sans-serif")
      .attr("dx", 3)
    let x;

    // Continuous
    if (color.interpolate) {
      const n = Math.min(color.domain().length, color.range().length);

      x = color.copy().rangeRound(d3.quantize(d3.interpolate(height - marginBottom, marginTop), n));

      svg.append("image")
          .attr("x", offsetLeft + marginLeft)
          .attr("y", offsetTop + marginTop)
          .attr("width", width - marginLeft - marginRight)
          .attr("height", height - marginTop - marginBottom)
          .attr("preserveAspectRatio", "none")
          .attr("xlink:href", ramp(color.copy().domain(d3.quantize(d3.interpolate(0, 1), n))).toDataURL());
    }

    // Sequential
    else if (color.interpolator) {
      x = Object.assign(color.copy()
          .interpolator(d3.interpolateRound(height - marginBottom, marginTop)),
          {range() { return [height - marginBottom, marginTop]; }});

      svg.append("image")
          .attr("x", offsetLeft + marginLeft)
          .attr("y", offsetTop + marginTop)
          .attr("width", width - marginLeft - marginRight)
          .attr("height", height - marginTop - marginBottom)
          .attr("preserveAspectRatio", "none")
          .attr("xlink:href", ramp(color.interpolator()).toDataURL());

      // scaleSequentialQuantile doesn’t implement ticks or tickFormat.
      if (!x.ticks) {
        if (tickValues === undefined) {
          const n = Math.round(ticks + 1);
          tickValues = d3.range(n).map(i => d3.quantile(color.domain(), i / (n - 1)));
        }
        if (typeof tickFormat !== "function") {
          tickFormat = d3.format(tickFormat === undefined ? ",f" : tickFormat);
        }
      }
    }

    // Threshold
    else if (color.invertExtent) {
      const thresholds
          = color.thresholds ? color.thresholds() // scaleQuantize
          : color.quantiles ? color.quantiles() // scaleQuantile
          : color.domain(); // scaleThreshold

      const thresholdFormat
          = tickFormat === undefined ? d => d
          : typeof tickFormat === "string" ? d3.format(tickFormat)
          : tickFormat;

      x = d3.scaleLinear()
          .domain([-1, color.range().length - 1])
          .rangeRound([height - marginBottom, marginTop]);

      svg.append("g")
        .selectAll("rect")
        .data(color.range())
        .join("rect")
          .attr("y", (d, i) => x(i))
          .attr("x", offsetLeft + marginLeft)
          .attr("height", (d, i) => x(i - 1) - x(i))
          .attr("width", width - marginRight - marginLeft)
          .attr("fill", d => d);

      tickValues = d3.range(thresholds.length);
      tickFormat = i => thresholdFormat(thresholds[i], i);
    }

    // Ordinal
    else {
      x = d3.scaleBand()
          .domain(color.domain())
          .rangeRound([height - marginBottom, marginTop]);

      svg.append("g")
        .selectAll("rect")
        .data(color.domain())
        .join("rect")
          .attr("y", x)
          .attr("x", offsetLeft + marginLeft)
          .attr("height", Math.max(0, x.bandwidth() - 1))
          .attr("width", width - marginLeft - marginRight)
          .attr("fill", color);

      tickAdjust = () => {};
    }

    svg.append("g")
        .attr("transform", `translate(${offsetLeft + width - marginRight},${offsetTop})`)
        .call(d3.axisRight(x)
          .ticks(ticks, typeof tickFormat === "string" ? tickFormat : undefined)
          .tickFormat(typeof tickFormat === "function" ? tickFormat : undefined)
          .tickSize(tickSize)
          .tickValues(tickValues))
        .call(tickAdjust)
        .call(g => g.select(".domain").remove())
        .call(g => g.append("text")
          .attr("x", marginLeft - width + marginRight)
          .attr("y", 0)
          .attr("fill", "currentColor")
          .attr("text-anchor", "start")
          .attr("font-weight", "bold")
          .attr("class", "title")
          .attr("font-family", "sans-serif")
          .attr("font-size", fontSize)
          .attr("font-weight", "bold")
          .text(title));

    return svg.node();
  }

	// Gene information

	function getGeneInfo(gene) {
    if (!state.query_geneinfo) {
      state.query_geneinfo = true
      console.log(`getGeneInfo("${gene}")`)
      mygene_query(gene, function(data) {
        console.log(`getGeneInfo("${gene}") > mygene_query`)
         g_blob = data
        if (Object.keys(data).indexOf("hits") !== -1) {
          if (data.hits[0]) {
            state.query_geneinfo = false
            var id = data.hits[0]._id;
            mygene_gene(id, fill_geneinfo);
          }
        }
      });
    }
  }

	function search_pubmed(query) {
		var url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?";
		$.get(url, {
			usehistory: "y",
			db: "pubmed",
			term: query
		})
		.done(function(data) {
			console.log("esearch done");
			var xml = $(data);
			var web = xml.find("WebEnv").html();
			var key = xml.find("QueryKey").html();
			var pubmed_ids = xml.find("IdList").children().map(function(i, item) {
				return item.innerHTML;
			});
			add_abstracts(query, pubmed_ids);
		});
	}

	function add_abstracts(query, pubmed_ids) {
		var url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?";
		$.get(url, {
			db: "pubmed",
			rettype: "abstract",
			id: Array.prototype.join.call(pubmed_ids)
		})
		.done(function(data) {
			console.log("efetch done");
			// Example: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=17284678,9997&rettype=abstract
			// x = data;
			var titles = $(data).find("ArticleTitle").map(function(i, item) {
				return {
					pubmed: pubmed_ids[i],
					title: item.innerHTML
				};
			// }).toArray().sort(dynamicSort('title')).map(function(item) {
			}).toArray().map(function(item) {
				return '<li><a href="https://www.ncbi.nlm.nih.gov/pubmed/' +
					item.pubmed + '">' + item.title + '</a></li>';
			}).join('');
			$('#pubmed-abstracts').innerHTML = '';
			$('<ul>' + titles + '</ul>').prependTo("#pubmed-abstracts");
		});
	}

	function mygene_query(gene, callback) {
		var request = new XMLHttpRequest();
		request.open(
			'GET', 'https://mygene.info/v3/query?q=symbol:' + gene.toLowerCase() + '&species=human', true
		);
		request.onload = function() {
			if (request.status >= 200 && request.status < 400) {
				// Success!
				var data = JSON.parse(request.responseText);
				callback(data);
			} else {
				// We reached our target server, but it returned an error
				console.log('error ' + request.responseText);
			}
		};
		request.onerror = function(data) {
			console.log('error ' + data);
		};
		request.send();
	}

	function mygene_gene(id, callback) {
		var request = new XMLHttpRequest();
		request.open(
			'GET',
			'https://mygene.info/v3/gene/' + id + '?email=kslowikowski@gmail.com'
		);
		request.onload = function() {
			if (request.status >= 200 && request.status < 400) {
				// Success!
				var data = JSON.parse(request.responseText);
				callback(data);
			} else {
				// We reached our target server, but it returned an error
				console.log('error ' + request.responseText);
			}
		};
		request.onerror = function(data) {
			console.log('error ' + data);
		};
		request.send();
	}

	function fill_geneinfo(data) {
		console.log(`fill_geneinfo()`)
		console.log(data)
		data.pubmed_term = 'immunotherapy'
		var geneinfo = document.getElementById('geneinfo')
	//  // var el = document.createElement("div");
		data.aliases = Array.isArray(data.alias) ? data.alias.join(", ") : data.alias
		//data.generif_list = data.generif.sort(dynamicSort('text')).map(function(d) {
		data.generif_list = data.generif.map(function(d) {
			return '<a href="https://www.ncbi.nlm.nih.gov/pubmed/' + d.pubmed + '">' +
				d.text + '</a>'
		}).map(function(d) {
			return '<li>' + d + '</li>'
		}).join('\n')
		// var form = '<h4><i>{symbol}</i></h4>' +
		// '<p>{name}</p>' +
		// '<p><b>Type:</b> {type_of_gene} <b>Aliases:</b> {aliases}</p>' +
		// '<p><b>Summary:</b> {summary}</p>' +
		// // '<p><b>iHOP:</b> <a target="_blank" href="http://www.ihop-net.org/UniPub/iHOP/index.html?field=all&search={symbol}&organism_id=1">link</a></p>' +
		// '<p><b>GeneRIF:</b><div id="generif"><ul>{generif_list}</ul></div></p>' +
		// '<p><b>Pubmed:</b> "{pubmed_term}" AND {symbol}' +
		// '<div id="pubmed-abstracts"></div></p>'
		// geneinfo.innerHTML = format(form, data)
		geneinfo.innerHTML = `<h4><i>${data.symbol}</i></h4>
		<p>${data.name}</p>
		<p><b>Type:</b> ${data.type_of_gene} <b>Aliases:</b> ${data.aliases}</p>
		<p><b>Summary:</b> ${data.summary}</p>
		<p><b>GeneRIF:</b> (${data.generif.length} results)<div id="generif"><ul>${data.generif_list}</ul></div></p>
		<p><b>Pubmed:</b> "${data.pubmed_term}" AND ${data.symbol}
		<div id="pubmed-abstracts"></div></p>`
		// search_pubmed(`"${data.pubmed_term}" AND ` + data.symbol)
	}

	/****************************************************************************
	 * Pure functions.
	 */

	/**
	 * Sort objects by a single property shared by all of the objects.
	 *
	 * @param {String} property - The property present in each object.
	 *
	 * @returns {Function} - A function to be use with Array.sort()
	 */
	function dynamicSort(property) {
			var sortOrder = 1;
			if (property[0] === "-") {
					sortOrder = -1;
					property = property.substr(1);
			}
			return function (a, b) {
					var result = (a[property] < b[property]) 
					? -1
					: (a[property] > b[property])
						? 1
						: 0;
							return result * sortOrder;
			};
	}

	/**
	 * Insert commas into a number so triplets are separated. Returns a string.
	 *
	 * @example
	 * // returns "1,000"
	 * numberWithCommas(1000);
	 *
	 * @param {Number} x - A number.
	 *
	 * @returns {String} - The number coverted to a string with commas.
	 */
	function numberWithCommas(x) {
		return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
	}

	/**
	 * Format a template string.
	 *
	 * @example
	 * // returns "Hi Kamil"
	 * format("Hi {name}", {name: "Kamil"});
	 *
	 * @param {String} form - The template string.
	 * @param {Object} datum - An object. Its properties are used in the template.
	 *
	 * @returns {String} - The formatted template with values filled in.
	 */
	function format(form, datum) {
		return form.replace(/{([^}]+)}/g, function(match, key) { 
			return typeof datum[key] != 'undefined' ? datum[key] : '';
		});
	}

	/**
	 * Convert an object with genomic coordinates into a region string.
	 *
	 * @example
	 * // returns "chr1:1,000-10,000,000"
	 * posToRegion({chr: "chr1", start: 1000, end: 10000000});
	 *
	 * @param {Object} pos - An object with properties: chr, start, end
	 * 
	 * @returns {String} - A string representation of the region.
	 */
	function posToRegion(pos) {
		if (pos['chr'] && pos['start'] && pos['end']) {
			return pos['chr'] + ':' +
				numberWithCommas(pos['start']) + '-' +
				numberWithCommas(pos['end']);
		}
		return '';
	}

	/**
	 * Get a nested property inside an object if it exists.
	 *
	 * @example
	 * // returns "deep secret"
	 * getNestedKey({html: body: { heart: "deep secret" } }, "html.body.heart");
	 *
	 * @example
	 * // returns null
	 * getNestedKey({html: body: { heart: "deep secret" } }, "html.body.brain");
	 *
	 * @param {Object} item - The object with nested properties.
	 * @param {String} keys - A period-delimited string of properties.
	 *
	 * @returns {Object} - The value stored inside the given object.
	 *
	 */
	function getNestedKey(item, keys) {
		var result = item;
		keys = keys.split(".");
		for (var i = 0; i < keys.length; i++) {
		if (result[keys[i]]) {
			result = result[keys[i]];
		} else {
			return null;
		}
		}
		return result;
	}

	/**
	 * Encode HTML entitites in a string.
	 *
	 * @example
	 * // returns "&#x3C;strong&#x3E;"
	 * encodeStr("<strong>");
	 *
	 * @param {String} raw - A string with HTML entities.
	 *
	 * @returns {String} - A string with HTML entities encoded.
	 */
	function encodeStr(raw) {
		return raw.replace(/[\u00A0-\u9999<>\&]/gim, function(i) {
		return '&#' + i.charCodeAt(0) + ';';
		});
	}

	/****************************************************************************
	 * Links to other websites.
	 */

	/**
	 * Create a link to the NCBI Entrez Gene database.
	 *
	 * @example
	 * // returns '<strong>Entrez:</strong> <a target="_blank" href="http://www.ncbi.nlm.nih.gov/gene/123">123</a>'
	 * entrezGeneLink({entrezgene: 123});
	 *
	 * @param {Object} datum - An object with a property: entrezgene
	 *
	 * @returns {String} - An HTML link to a gene at the NCBI Entrez Gene database.
	 */
	function entrezGeneLink(datum) {
		var form = '<strong>Entrez:</strong> ' +
			'<a target="_blank" href="http://www.ncbi.nlm.nih.gov/gene/{entrezgene}">' +
			'{entrezgene}</a>';
		if (datum['entrezgene']) {
			return format(form, datum);
		}
		return '';
	}

	/**
	 * Create a link to the HGNC database.
	 *
	 * @example
	 * // returns '<strong>HGNC:</strong> <a target="_blank" href="http://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:123">123</a>'
	 * hgncGeneLink({HGNC: 123});
	 *
	 * @param {Object} datum - An object with a property: HGNC
	 *
	 * @returns {String} - An HTML link to a gene at the HGNC database.
	 */
	function hgncGeneLink(datum) {
		var form = '<a target="_blank"' +
			' href="http://www.genenames.org/cgi-bin/gene_symbol_report' +
			'?hgnc_id=HGNC:{HGNC}">HGNC</a>';
		if (datum['HGNC']) {
			return format(form, datum);
		}
		return '';
	}

	/**
	 * Create an HTML list of pathways described in WikiPathways.
	 *
	 * @param {Object} datum - An object with a property: pathway.wikipathways
	 *
	 * @returns {String} - An HTML list of pathways.
	 */
	function wikipathwayLinks(datum) {
		var links = '<strong>Wikipathways:</strong>';
		var pathways = datum['pathway.wikipathways'];
		var url = 'http://www.wikipathways.org/index.php/Pathway:{id}';
		var form = '<li><a target="_blank" href="' + url + '">{name}</a></li><br>';
		if (pathways) {
			links += '<ul>';
			for (var i = 0; i < pathways.length; i++) {
				links += format(form, pathways[i]);
			}
			links += '</ul>';
		}
		return links;
	}

	/**
	 * Create an HTML list of the sentences in a gene's summary.
	 *
	 * @example
	 * //returns "<ul><li>One.</li><li>Two.</li></ul>"
	 * summaryList({summary: "One. Two."});
	 *
	 * @param {Object} datum - An object with property: datum
	 *
	 * @returns {String} - An HTML list.
	 */
	function summaryList(datum) {
		var summary = datum['summary'];
		var result = '';
		if (summary) {
			var parts = summary.split(".");
			result += '<ul>';
			for (var i = 0; i < parts.length - 1; i++) {
				result +=  '<li>' + parts[i] + '.</li><br>';
			}
			result += '</ul>';
		}
		return result;
	}

	/**
	 * Create a link to the UCSC genome browser.
	 *
	 * @param {Object} pos - An object with properties: chr, start, end
	 *
	 * @returns {String} - The url to the UCSC genome browser.
	 */
	function ucscRegionLink(pos) {
		var form = '<a target="_blank"' +
			' href="https://genome.ucsc.edu/cgi-bin/hgTracks' +
			'?db=hg19&position=chr{chr}%3A{start}-{end}">' +
			'UCSC</a>';
		if (pos['chr'] && pos['start'] && pos['end']) {
			return format(form, pos);
		}
		return '';
	}


  // Only export these functions
  return {
    "main": main,
    "onClusterNameClick": onClusterNameClick,
    "onGeneGroupbyClick": onGeneGroupbyClick,
    "onMarkerGeneClick": onMarkerGeneClick,
    "renderer": renderer,
    "state": state
  }


}();
