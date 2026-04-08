// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, true, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "Ofav_080", "Ofav_080", "Emerald", "1", "1", "Eo1-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:38", "1", "0.513999652", "0.327", "0.582", "10", "1", "0", "22.569637", "8.066104" ], [ "2", "Ofav_081", "Ofav_081", "Emerald", "3", "1", "Eo3-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:40", "3", "0.332052499", "0.434", "0.737", "18", "0.93", "0.07", "22.517529", "8.065329" ], [ "3", "Ofav_082", "Ofav_082", "Emerald", "4", "1", "Eo4-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:42", "4", "0.437211206", "0.294", "0.507", "31", "1", "0", "22.523381", "8.067308" ], [ "4", "Ofav_083", "Ofav_083", "Emerald", "5", "1", "Eo5-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:43", "5", "0.548562048", "0.358", "0.606", "10", "1", "0", "22.525446", "8.067308" ], [ "5", "Ofav_084", "Ofav_084", "Rainbow", "3", "1", "Ro3-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:46", "8", "0.683454293", "0.323", "0.551", "17", "1", "0", "22.530609", "8.064591" ], [ "6", "Ofav_085", "Ofav_085", "Rainbow", "4", "1", "Ro4-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:47", "9", "0.420509784", "0.388", "0.636", "17", "0.9", "0.11", "22.537625", "8.065505" ], [ "7", "Ofav_086", "Ofav_086", "Rainbow", "5", "1", "Ro5-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:48", "10", "0.42361281", "0.397", "0.657", "32", "0", "1", "22.545034", "8.067564" ], [ "8", "Ofav_087", "Ofav_087", "Star", "3", "1", "So3-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:51", "13", "0.265474824", "0.43", "0.72", "25", "0", "1", "22.555262", "8.065163" ], [ "9", "Ofav_088", "Ofav_088", "Star", "4", "1", "So4-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:52", "14", "0.269329097", "0.437", "0.739", "31", "0", "1", "22.571604", "8.065336" ], [ "10", "Ofav_089", "Ofav_089", "Star", "5", "1", "So5-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:53", "15", "0.209789356", "0.486", "0.758", "31", "0", "1", "22.570833", "8.066034" ], [ "11", "Ofav_090", "Ofav_090", "MacN", "1", "1", "Mo1-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:54", "16", "-0.132082631", "0.402", "0.722", "28", "0", "1", "22.592339", "8.067567" ], [ "12", "Ofav_091", "Ofav_091", "MacN", "4", "1", "Mo4-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:56", "19", "-0.105226315", "0.429", "0.737", "31", "0", "1", "22.598748", "8.066505" ], [ "13", "Ofav_092", "Ofav_092", "MacN", "5", "1", "Mo5-1", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:57", "20", "0.288952665", "0.485", "0.808", "34", "0", "1", "22.592601", "8.065325" ], [ "14", "Ofav_093", "Ofav_093", "MacN", "6", "1", "CR043", "1", "controlpH", "controltemp", "controlpH_controltemp", "CC", "02:58", "21", "0.038396602", "0.48", "0.782", "25", "0", "1", "22.603895", "8.064723" ], [ "15", "Ofav_094", "Ofav_094", "Emerald", "1", "2", "Eo1-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:01", "1", "0.393152124", "0.342", "0.631", "12", "1", "0", "33.004931", "7.745075" ], [ "16", "Ofav_095", "Ofav_095", "Emerald", "2", "2", "Eo2-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:02", "2", "0.199876607", "0.452", "0.755", "24", "0.3", "0.7", "32.991769", "7.743074" ], [ "17", "Ofav_096", "Ofav_096", "Emerald", "3", "2", "Eo3-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:03", "3", "0.211774423", "0.447", "0.729", "23", "0.67", "0.33", "32.974558", "7.748701" ], [ "18", "Ofav_097", "Ofav_097", "Emerald", "4", "2", "Eo4-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:05", "4", "0.222988735", "0.314", "0.535", "32", "1", "0", "32.95538", "7.7615" ], [ "19", "Ofav_098", "Ofav_098", "Emerald", "5", "2", "Eo5-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:06", "5", "0.403904178", "0.366", "0.618", "11", "0.92", "0.08", "32.954085", "7.773905" ], [ "20", "Ofav_099", "Ofav_099", "Rainbow", "1", "2", "Ro1-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:09", "6", "-0.028996722", "0.377", "0.658", "13", "0.63", "0.38", "32.96156", "7.836482" ], [ "21", "Ofav_100", "Ofav_100", "Rainbow", "2", "2", "Ro2-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:10", "7", "0.056504208", "0.359", "0.653", "14", "0", "1", "32.965281", "7.852662" ], [ "22", "Ofav_101", "Ofav_101", "Rainbow", "3", "2", "Ro3-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:10", "8", "0.265040369", "0.325", "0.553", "11", "0", "0", "32.965281", "7.852662" ], [ "23", "Ofav_102", "Ofav_102", "Rainbow", "4", "2", "Ro4-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:08", "9", "0.327051796", "0.401", "0.652", "11", "1", "0", "32.954446", "7.823464" ], [ "24", "Ofav_103", "Ofav_103", "Rainbow", "5", "2", "Ro5-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:07", "10", "0.334279926", "0.365", "0.649", "25", "0", "1", "32.953872", "7.785509" ], [ "25", "Ofav_104", "Ofav_104", "Star", "1", "2", "So1-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:13", "11", "0.074842059", "0.414", "0.721", "16", "0", "1", "32.999112", "7.957457" ], [ "26", "Ofav_105", "Ofav_105", "Star", "2", "2", "So2-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:15", "12", "-0.026005504", "0.507", "0.827", "18", "0", "1", "33.017225", "8.004322" ], [ "27", "Ofav_106", "Ofav_106", "Star", "3", "2", "So3-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:16", "13", "0.29230477", "0.464", "0.744", "14", "0", "1", "33.022289", "8.026373" ], [ "28", "Ofav_107", "Ofav_107", "Star", "4", "2", "So4-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:17", "14", "0.141618714", "0.467", "0.804", "21", "0", "1", "33.019454", "8.028577" ], [ "29", "Ofav_108", "Ofav_108", "Star", "5", "2", "So5-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:19", "15", "0.424737112", "0.421", "0.683", "30", "0", "1", "32.988458", "7.997266" ], [ "30", "Ofav_109", "Ofav_109", "MacN", "1", "2", "Mo1-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:22", "16", "-0.055451417", "0.429", "0.773", "18", "0", "1", "32.945152", "7.890439" ], [ "31", "Ofav_110", "Ofav_110", "MacN", "2", "2", "Mo2-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:25", "17", "-0.080435075", "0.399", "0.759", "11", "0", "1", "32.950774", "7.851446" ], [ "32", "Ofav_111", "Ofav_111", "MacN", "3", "2", "Mo3-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:24", "18", "0.22596648", "0.465", "0.75", "35", "0", "1", "32.943398", "7.859227" ], [ "33", "Ofav_112", "Ofav_112", "MacN", "4", "2", "Mo4-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:28", "19", "0.132260971", "0.439", "0.78", "28", "0", "1", "32.972689", "7.845783" ], [ "34", "Ofav_113", "Ofav_113", "MacN", "5", "2", "Mo5-2", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:24", "20", "0.259780238", "0.476", "0.812", "25", "0", "1", "32.943398", "7.859227" ], [ "35", "Ofav_114", "Ofav_114", "MacN", "6", "2", "CR046", "6", "controlpH", "hightemp", "controlpH_hightemp", "CH", "03:22", "21", "0.298560037", "0.492", "0.812", "35", "0", "1", "32.945152", "7.890439" ], [ "36", "Ofav_115", "Ofav_115", "Emerald", "1", "3", "Eo1-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:32", "1", "0.229527512", "0.351", "0.601", "NA", "1", "0", "22.505449", "7.815671" ], [ "37", "Ofav_116", "Ofav_116", "Emerald", "3", "3", "Eo3-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:34", "3", "0.233784935", "0.317", "0.547", "NA", "0.97", "0.03", "22.501449", "7.814456" ], [ "38", "Ofav_117", "Ofav_117", "Emerald", "4", "3", "Eo4-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:35", "4", "0.391298251", "0.295", "0.516", "17", "1", "0", "22.506071", "7.815627" ], [ "39", "Ofav_118", "Ofav_118", "Emerald", "5", "3", "Eo5-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:36", "5", "0.562448949", "0.336", "0.574", "10", "1", "0", "22.518267", "7.822422" ], [ "40", "Ofav_119", "Ofav_119", "Rainbow", "2", "3", "Ro2-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:40", "7", "0.132500708", "0.322", "0.593", "NA", "0", "1", "22.54782", "7.819187" ], [ "41", "Ofav_120", "Ofav_120", "Rainbow", "4", "3", "Ro4-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:41", "9", "0.181406898", "0.409", "0.683", "NA", "0.73", "0.27", "22.559458", "7.818991" ], [ "42", "Ofav_121", "Ofav_121", "Rainbow", "5", "3", "Ro5-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:42", "10", "0.274811332", "0.367", "0.637", "25", "0", "1", "22.564834", "7.818991" ], [ "43", "Ofav_122", "Ofav_122", "Star", "1", "3", "So1-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:45", "11", "0.009955487", "0.367", "0.637", "NA", "0", "1", "22.577505", "7.815602" ], [ "44", "Ofav_123", "Ofav_123", "Star", "2", "3", "So2-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:46", "12", "0.096023941", "0.405", "0.66", "25", "0", "1", "22.592847", "7.820235" ], [ "45", "Ofav_124", "Ofav_124", "Star", "3", "3", "So3-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:47", "13", "0.333742882", "0.417", "0.703", "25", "0", "1", "22.597813", "7.820235" ], [ "46", "Ofav_125", "Ofav_125", "Star", "4", "3", "So4-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:49", "14", "0.079499847", "0.423", "0.703", "27", "0", "1", "22.604911", "7.822518" ], [ "47", "Ofav_126", "Ofav_126", "Star", "5", "3", "So5-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:50", "15", "0.350097769", "0.414", "0.655", "27", "0", "1", "22.59724", "7.822643" ], [ "48", "Ofav_127", "Ofav_127", "MacN", "1", "3", "Mo1-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:47", "16", "-0.183405798", "0.376", "0.695", "27", "0", "1", "22.597813", "7.820235" ], [ "49", "Ofav_128", "Ofav_128", "MacN", "4", "3", "Mo4-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:48", "19", "-0.021731628", "0.439", "0.733", "27", "0", "1", "22.606239", "7.821077" ], [ "50", "Ofav_129", "Ofav_129", "MacN", "5", "3", "Mo5-3", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:51", "20", "0.195532221", "0.457", "0.746", "NA", "0", "1", "22.55577", "7.823256" ], [ "51", "Ofav_130", "Ofav_130", "MacN", "6", "3", "CR047", "4", "lowpH", "controltemp", "lowpH_controltemp", "LC", "03:52", "21", "0.043069287", "0.421", "0.697", "27", "0", "1", "22.521676", "7.823712" ], [ "52", "Ofav_131", "Ofav_131", "Emerald", "2", "4", "Eo2-4", "8", "lowpH", "hightemp", "lowpH_hightemp", "LH", "03:58", "2", "0.420727629", "0.425", "0.725", "14", "1", "0", "32.995916", "7.821764" ], [ "53", "Ofav_132", "Ofav_132", "Emerald", "4", "4", "Eo4-4", "8", "lowpH", "hightemp", "lowpH_hightemp", "LH", "04:00", "4", "0.418455126", "0.421", "0.687", "21", "1", "0", "32.993506", "7.823659" ], [ "54", "Ofav_133", "Ofav_133", "Emerald", "5", "4", "Eo5-4", "8", "lowpH", "hightemp", "lowpH_hightemp", "LH", "04:01", "5", "0.668131785", "0.399", "0.669", "7", "1", "0", "32.988245", "7.823673" ], [ "55", "Ofav_153", "Ofav_153", "Emerald", "4", "3", "Eo4-3", "4", "lowpH", "hightemp", "lowpH_hightemp", "LH", "03:05", "4", "0.391298251", "0.22", "0.385", "17", "1", "0", "32.998965", "7.815966" ], [ "56", "Ofav_154", "Ofav_154", "Emerald", "5", "3", "Eo5-3", "4", "lowpH", "hightemp", "lowpH_hightemp", "LH", "03:07", "5", "0.562448949", "0.329", "0.562", "10", "0", "0", "33.024306", "7.816174" ], [ "57", "Ofav_134", "Ofav_134", "Rainbow", "5", "4", "Ro5-4", "8", "lowpH", "hightemp", "lowpH_hightemp", "LH", "04:02", "10", "0.565510598", "0.388", "0.643", "25", "0", "1", "32.994162", "7.822569" ], [ "58", "Ofav_155", "Ofav_155", "Rainbow", "5", "3", "Ro5-3", "4", "lowpH", "hightemp", "lowpH_hightemp", "LH", "03:09", "10", "0.274811332", "0.368", "0.639", "25", "0", "1", "33.033698", "7.817851" ], [ "59", "Ofav_135", "Ofav_135", "Star", "5", "4", "So5-4", "8", "lowpH", "hightemp", "lowpH_hightemp", "LH", "04:03", "15", "0.324732873", "0.47", "0.731", "25", "0", "1", "32.995047", "7.822829" ], [ "60", "Ofav_156", "Ofav_156", "Star", "2", "3", "So2-3", "4", "lowpH", "hightemp", "lowpH_hightemp", "LH", "03:16", "12", "0.096023941", "0.394", "0.642", "25", "0", "1", "32.957462", "7.823444" ], [ "61", "Ofav_157", "Ofav_157", "Star", "3", "3", "So3-3", "4", "lowpH", "hightemp", "lowpH_hightemp", "LH", "03:17", "13", "0.333742882", "0.396", "0.668", "25", "0", "1", "32.94966", "7.823404" ], [ "62", "Ofav_158", "Ofav_158", "Star", "4", "3", "So4-3", "4", "lowpH", "hightemp", "lowpH_hightemp", "LH", "03:18", "14", "0.079499847", "0.415", "0.689", "27", "0", "1", "32.937333", "7.823404" ], [ "63", "Ofav_159", "Ofav_159", "Star", "5", "3", "So5-3", "4", "lowpH", "hightemp", "lowpH_hightemp", "LH", "03:19", "15", "0.350097769", "0.383", "0.606", "27", "0", "1", "32.933236", "7.824068" ], [ "64", "Ofav_136", "Ofav_136", "MacN", "4", "4", "Mo4-4", "8", "lowpH", "hightemp", "lowpH_hightemp", "LH", "04:05", "19", "0.148192161", "0.497", "0.844", "24", "0", "1", "33.027059", "7.823024" ], [ "65", "Ofav_137", "Ofav_137", "MacN", "5", "4", "Mo5-4", "8", "lowpH", "hightemp", "lowpH_hightemp", "LH", "04:05", "20", "0.063717664", "0.491", "0.792", "25", "0", "1", "33.027059", "7.823024" ], [ "66", "Ofav_138", "Ofav_138", "MacN", "6", "4", "CR049", "8", "lowpH", "hightemp", "lowpH_hightemp", "LH", "04:03", "21", "1.140602553", "0.482", "0.782", "25", "0", "1", "32.995047", "7.822829" ], [ "67", "Ofav_160", "Ofav_160", "MacN", "1", "3", "Mo1-3", "4", "lowpH", "hightemp", "lowpH_hightemp", "LH", "03:24", "16", "-0.183405798", "0.392", "0.725", "27", "0", "1", "32.961953", "7.822263" ], [ "68", "Ofav_161", "Ofav_161", "MacN", "4", "3", "Mo4-3", "4", "lowpH", "hightemp", "lowpH_hightemp", "LH", "03:25", "19", "-0.021731628", "0.417", "0.696", "27", "0", "1", "32.969001", "7.822056" ], [ "69", "Ofav_162", "Ofav_162", "MacN", "5", "3", "Mo5-3", "4", "lowpH", "hightemp", "lowpH_hightemp", "LH", "03:26", "20", "0.195532221", "0.419", "0.684", "NA", "0", "1", "32.989687", "7.822228" ], [ "70", "Ofav_163", "Ofav_163", "MacN", "6", "3", "CR047", "4", "lowpH", "hightemp", "lowpH_hightemp", "LH", "03:27", "21", "0.043069287", "0.429", "0.71", "27", "0", "1", "32.993818", "7.822276" ] ];
var svgObjectNames   = [ "pca", "dens" ];

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
