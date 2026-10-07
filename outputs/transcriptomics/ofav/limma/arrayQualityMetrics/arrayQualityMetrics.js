// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = ["false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","true","false","false","false","false","true","false","false","false","false","false","false","false"];
var arrayMetadata    = [{"array":"1","sampleNames":"Ofav_080","id":"Ofav_080","site":"Emerald","geno":"1","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:38","genotype":"1","propB":"1","propD":"0","temp_actual":"22.569637","ph_actual":"8.066104","habitat":"reef"},{"array":"2","sampleNames":"Ofav_081","id":"Ofav_081","site":"Emerald","geno":"3","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:40","genotype":"3","propB":"0.93","propD":"0.07","temp_actual":"22.517529","ph_actual":"8.065329","habitat":"reef"},{"array":"3","sampleNames":"Ofav_082","id":"Ofav_082","site":"Emerald","geno":"4","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:42","genotype":"4","propB":"1","propD":"0","temp_actual":"22.523381","ph_actual":"8.067308","habitat":"reef"},{"array":"4","sampleNames":"Ofav_083","id":"Ofav_083","site":"Emerald","geno":"5","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:43","genotype":"5","propB":"1","propD":"0","temp_actual":"22.525446","ph_actual":"8.067308","habitat":"reef"},{"array":"5","sampleNames":"Ofav_084","id":"Ofav_084","site":"Rainbow","geno":"3","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:46","genotype":"8","propB":"1","propD":"0","temp_actual":"22.530609","ph_actual":"8.064591","habitat":"reef"},{"array":"6","sampleNames":"Ofav_085","id":"Ofav_085","site":"Rainbow","geno":"4","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:47","genotype":"9","propB":"0.9","propD":"0.11","temp_actual":"22.537625","ph_actual":"8.065505","habitat":"reef"},{"array":"7","sampleNames":"Ofav_086","id":"Ofav_086","site":"Rainbow","geno":"5","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:48","genotype":"10","propB":"0","propD":"1","temp_actual":"22.545034","ph_actual":"8.067564","habitat":"reef"},{"array":"8","sampleNames":"Ofav_087","id":"Ofav_087","site":"Star","geno":"3","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:51","genotype":"13","propB":"0","propD":"1","temp_actual":"22.555262","ph_actual":"8.065163","habitat":"urban"},{"array":"9","sampleNames":"Ofav_088","id":"Ofav_088","site":"Star","geno":"4","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:52","genotype":"14","propB":"0","propD":"1","temp_actual":"22.571604","ph_actual":"8.065336","habitat":"urban"},{"array":"10","sampleNames":"Ofav_089","id":"Ofav_089","site":"Star","geno":"5","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:53","genotype":"15","propB":"0","propD":"1","temp_actual":"22.570833","ph_actual":"8.066034","habitat":"urban"},{"array":"11","sampleNames":"Ofav_090","id":"Ofav_090","site":"MacN","geno":"1","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:54","genotype":"16","propB":"0","propD":"1","temp_actual":"22.592339","ph_actual":"8.067567","habitat":"urban"},{"array":"12","sampleNames":"Ofav_091","id":"Ofav_091","site":"MacN","geno":"4","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:56","genotype":"19","propB":"0","propD":"1","temp_actual":"22.598748","ph_actual":"8.066505","habitat":"urban"},{"array":"13","sampleNames":"Ofav_092","id":"Ofav_092","site":"MacN","geno":"5","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:57","genotype":"20","propB":"0","propD":"1","temp_actual":"22.592601","ph_actual":"8.065325","habitat":"urban"},{"array":"14","sampleNames":"Ofav_093","id":"Ofav_093","site":"MacN","geno":"6","replicate":"1","tank":"1","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"02:58","genotype":"21","propB":"0","propD":"1","temp_actual":"22.603895","ph_actual":"8.064723","habitat":"urban"},{"array":"15","sampleNames":"Ofav_094","id":"Ofav_094","site":"Emerald","geno":"1","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:01","genotype":"1","propB":"1","propD":"0","temp_actual":"33.004931","ph_actual":"7.745075","habitat":"reef"},{"array":"16","sampleNames":"Ofav_095","id":"Ofav_095","site":"Emerald","geno":"2","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:02","genotype":"2","propB":"0.3","propD":"0.7","temp_actual":"32.991769","ph_actual":"7.743074","habitat":"reef"},{"array":"17","sampleNames":"Ofav_096","id":"Ofav_096","site":"Emerald","geno":"3","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:03","genotype":"3","propB":"0.67","propD":"0.33","temp_actual":"32.974558","ph_actual":"7.748701","habitat":"reef"},{"array":"18","sampleNames":"Ofav_097","id":"Ofav_097","site":"Emerald","geno":"4","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:05","genotype":"4","propB":"1","propD":"0","temp_actual":"32.95538","ph_actual":"7.7615","habitat":"reef"},{"array":"19","sampleNames":"Ofav_098","id":"Ofav_098","site":"Emerald","geno":"5","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:06","genotype":"5","propB":"0.92","propD":"0.08","temp_actual":"32.954085","ph_actual":"7.773905","habitat":"reef"},{"array":"20","sampleNames":"Ofav_099","id":"Ofav_099","site":"Rainbow","geno":"1","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:09","genotype":"6","propB":"0.63","propD":"0.38","temp_actual":"32.96156","ph_actual":"7.836482","habitat":"reef"},{"array":"21","sampleNames":"Ofav_100","id":"Ofav_100","site":"Rainbow","geno":"2","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:10","genotype":"7","propB":"0","propD":"1","temp_actual":"32.965281","ph_actual":"7.852662","habitat":"reef"},{"array":"22","sampleNames":"Ofav_101","id":"Ofav_101","site":"Rainbow","geno":"3","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:10","genotype":"8","propB":"0","propD":"0","temp_actual":"32.965281","ph_actual":"7.852662","habitat":"reef"},{"array":"23","sampleNames":"Ofav_102","id":"Ofav_102","site":"Rainbow","geno":"4","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:08","genotype":"9","propB":"1","propD":"0","temp_actual":"32.954446","ph_actual":"7.823464","habitat":"reef"},{"array":"24","sampleNames":"Ofav_103","id":"Ofav_103","site":"Rainbow","geno":"5","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:07","genotype":"10","propB":"0","propD":"1","temp_actual":"32.953872","ph_actual":"7.785509","habitat":"reef"},{"array":"25","sampleNames":"Ofav_104","id":"Ofav_104","site":"Star","geno":"1","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:13","genotype":"11","propB":"0","propD":"1","temp_actual":"32.999112","ph_actual":"7.957457","habitat":"urban"},{"array":"26","sampleNames":"Ofav_105","id":"Ofav_105","site":"Star","geno":"2","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:15","genotype":"12","propB":"0","propD":"1","temp_actual":"33.017225","ph_actual":"8.004322","habitat":"urban"},{"array":"27","sampleNames":"Ofav_106","id":"Ofav_106","site":"Star","geno":"3","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:16","genotype":"13","propB":"0","propD":"1","temp_actual":"33.022289","ph_actual":"8.026373","habitat":"urban"},{"array":"28","sampleNames":"Ofav_107","id":"Ofav_107","site":"Star","geno":"4","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:17","genotype":"14","propB":"0","propD":"1","temp_actual":"33.019454","ph_actual":"8.028577","habitat":"urban"},{"array":"29","sampleNames":"Ofav_108","id":"Ofav_108","site":"Star","geno":"5","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:19","genotype":"15","propB":"0","propD":"1","temp_actual":"32.988458","ph_actual":"7.997266","habitat":"urban"},{"array":"30","sampleNames":"Ofav_109","id":"Ofav_109","site":"MacN","geno":"1","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:22","genotype":"16","propB":"0","propD":"1","temp_actual":"32.945152","ph_actual":"7.890439","habitat":"urban"},{"array":"31","sampleNames":"Ofav_110","id":"Ofav_110","site":"MacN","geno":"2","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:25","genotype":"17","propB":"0","propD":"1","temp_actual":"32.950774","ph_actual":"7.851446","habitat":"urban"},{"array":"32","sampleNames":"Ofav_111","id":"Ofav_111","site":"MacN","geno":"3","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:24","genotype":"18","propB":"0","propD":"1","temp_actual":"32.943398","ph_actual":"7.859227","habitat":"urban"},{"array":"33","sampleNames":"Ofav_112","id":"Ofav_112","site":"MacN","geno":"4","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:28","genotype":"19","propB":"0","propD":"1","temp_actual":"32.972689","ph_actual":"7.845783","habitat":"urban"},{"array":"34","sampleNames":"Ofav_113","id":"Ofav_113","site":"MacN","geno":"5","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:24","genotype":"20","propB":"0","propD":"1","temp_actual":"32.943398","ph_actual":"7.859227","habitat":"urban"},{"array":"35","sampleNames":"Ofav_114","id":"Ofav_114","site":"MacN","geno":"6","replicate":"2","tank":"6","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"03:22","genotype":"21","propB":"0","propD":"1","temp_actual":"32.945152","ph_actual":"7.890439","habitat":"urban"},{"array":"36","sampleNames":"Ofav_117","id":"Ofav_117","site":"Emerald","geno":"4","replicate":"3","tank":"4","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"03:35","genotype":"4","propB":"1","propD":"0","temp_actual":"22.506071","ph_actual":"7.815627","habitat":"reef"},{"array":"37","sampleNames":"Ofav_118","id":"Ofav_118","site":"Emerald","geno":"5","replicate":"3","tank":"4","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"03:36","genotype":"5","propB":"1","propD":"0","temp_actual":"22.518267","ph_actual":"7.822422","habitat":"reef"},{"array":"38","sampleNames":"Ofav_121","id":"Ofav_121","site":"Rainbow","geno":"5","replicate":"3","tank":"4","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"03:42","genotype":"10","propB":"0","propD":"1","temp_actual":"22.564834","ph_actual":"7.818991","habitat":"reef"},{"array":"39","sampleNames":"Ofav_123","id":"Ofav_123","site":"Star","geno":"2","replicate":"3","tank":"4","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"03:46","genotype":"12","propB":"0","propD":"1","temp_actual":"22.592847","ph_actual":"7.820235","habitat":"urban"},{"array":"40","sampleNames":"Ofav_124","id":"Ofav_124","site":"Star","geno":"3","replicate":"3","tank":"4","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"03:47","genotype":"13","propB":"0","propD":"1","temp_actual":"22.597813","ph_actual":"7.820235","habitat":"urban"},{"array":"41","sampleNames":"Ofav_125","id":"Ofav_125","site":"Star","geno":"4","replicate":"3","tank":"4","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"03:49","genotype":"14","propB":"0","propD":"1","temp_actual":"22.604911","ph_actual":"7.822518","habitat":"urban"},{"array":"42","sampleNames":"Ofav_126","id":"Ofav_126","site":"Star","geno":"5","replicate":"3","tank":"4","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"03:50","genotype":"15","propB":"0","propD":"1","temp_actual":"22.59724","ph_actual":"7.822643","habitat":"urban"},{"array":"43","sampleNames":"Ofav_127","id":"Ofav_127","site":"MacN","geno":"1","replicate":"3","tank":"4","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"03:47","genotype":"16","propB":"0","propD":"1","temp_actual":"22.597813","ph_actual":"7.820235","habitat":"urban"},{"array":"44","sampleNames":"Ofav_128","id":"Ofav_128","site":"MacN","geno":"4","replicate":"3","tank":"4","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"03:48","genotype":"19","propB":"0","propD":"1","temp_actual":"22.606239","ph_actual":"7.821077","habitat":"urban"},{"array":"45","sampleNames":"Ofav_130","id":"Ofav_130","site":"MacN","geno":"6","replicate":"3","tank":"4","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"03:52","genotype":"21","propB":"0","propD":"1","temp_actual":"22.521676","ph_actual":"7.823712","habitat":"urban"},{"array":"46","sampleNames":"Ofav_131","id":"Ofav_131","site":"Emerald","geno":"2","replicate":"4","tank":"8","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"03:58","genotype":"2","propB":"1","propD":"0","temp_actual":"32.995916","ph_actual":"7.821764","habitat":"reef"},{"array":"47","sampleNames":"Ofav_132","id":"Ofav_132","site":"Emerald","geno":"4","replicate":"4","tank":"8","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"04:00","genotype":"4","propB":"1","propD":"0","temp_actual":"32.993506","ph_actual":"7.823659","habitat":"reef"},{"array":"48","sampleNames":"Ofav_133","id":"Ofav_133","site":"Emerald","geno":"5","replicate":"4","tank":"8","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"04:01","genotype":"5","propB":"1","propD":"0","temp_actual":"32.988245","ph_actual":"7.823673","habitat":"reef"},{"array":"49","sampleNames":"Ofav_134","id":"Ofav_134","site":"Rainbow","geno":"5","replicate":"4","tank":"8","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"04:02","genotype":"10","propB":"0","propD":"1","temp_actual":"32.994162","ph_actual":"7.822569","habitat":"reef"},{"array":"50","sampleNames":"Ofav_135","id":"Ofav_135","site":"Star","geno":"5","replicate":"4","tank":"8","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"04:03","genotype":"15","propB":"0","propD":"1","temp_actual":"32.995047","ph_actual":"7.822829","habitat":"urban"},{"array":"51","sampleNames":"Ofav_136","id":"Ofav_136","site":"MacN","geno":"4","replicate":"4","tank":"8","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"04:05","genotype":"19","propB":"0","propD":"1","temp_actual":"33.027059","ph_actual":"7.823024","habitat":"urban"},{"array":"52","sampleNames":"Ofav_137","id":"Ofav_137","site":"MacN","geno":"5","replicate":"4","tank":"8","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"04:05","genotype":"20","propB":"0","propD":"1","temp_actual":"33.027059","ph_actual":"7.823024","habitat":"urban"},{"array":"53","sampleNames":"Ofav_138","id":"Ofav_138","site":"MacN","geno":"6","replicate":"4","tank":"8","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"04:03","genotype":"21","propB":"0","propD":"1","temp_actual":"32.995047","ph_actual":"7.822829","habitat":"urban"}];
var svgObjectNames   = ["pca","dens"];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
