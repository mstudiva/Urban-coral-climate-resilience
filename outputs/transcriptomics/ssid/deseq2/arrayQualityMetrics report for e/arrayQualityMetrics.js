// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false ];
var arrayMetadata    = [ [ "1", "Ssid_001", "Ssid_001", "Emerald", "1", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "11:36", "22", "0.349", "0.759", "25", "0", "0.98", "22.576029", "8.041769" ], [ "2", "Ssid_002", "Ssid_002", "Emerald", "2", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "11:40", "23", "0.377", "0.818", "10", "0", "0", "22.606222", "8.040733" ], [ "3", "Ssid_003", "Ssid_003", "Emerald", "3", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "11:43", "24", "0.442", "0.811", "12", "0", "0.87", "22.567883", "8.044771" ], [ "4", "Ssid_004", "Ssid_004", "Emerald", "4", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "11:47", "25", "0.344", "0.709", "17", "0", "0.99", "22.490319", "8.045541" ], [ "5", "Ssid_005", "Ssid_005", "Emerald", "5", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "11:50", "26", "0.422", "0.752", "10", "0", "0.5", "22.487894", "8.045853" ], [ "6", "Ssid_006", "Ssid_006", "Rainbow", "1", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "11:55", "27", "0.364", "0.754", "25", "0", "1", "22.529544", "8.046972" ], [ "7", "Ssid_007", "Ssid_007", "Rainbow", "2", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "11:58", "28", "0.462", "0.885", "17", "1", "0", "22.551016", "8.047554" ], [ "8", "Ssid_008", "Ssid_008", "Rainbow", "3", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "12:01", "29", "0.321", "0.964", "5", "0", "0", "22.568932", "8.039805" ], [ "9", "Ssid_009", "Ssid_009", "Rainbow", "4", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "12:05", "30", "0.34", "0.653", "19", "0", "1", "22.602862", "8.046859" ], [ "10", "Ssid_010", "Ssid_010", "Rainbow", "5", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "12:08", "31", "0.363", "0.671", "25", "0", "1", "22.588241", "8.037466" ], [ "11", "Ssid_011", "Ssid_011", "Star", "1", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "12:07", "32", "0.34", "0.685", "19", "0", "1", "22.605321", "8.045583" ], [ "12", "Ssid_012", "Ssid_012", "Star", "2", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "12:13", "33", "0.395", "0.748", "12", "0", "1", "22.499794", "8.055674" ], [ "13", "Ssid_013", "Ssid_013", "Star", "3", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "12:17", "34", "0.371", "0.763", "19", "0", "1", "22.503498", "8.053427" ], [ "14", "Ssid_014", "Ssid_014", "Star", "4", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "12:18", "35", "0.302", "0.651", "11", "0", "1", "22.5113", "8.056571" ], [ "15", "Ssid_015", "Ssid_015", "Star", "5", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "12:19", "36", "0.428", "0.768", "27", "0", "1", "22.514956", "8.054615" ], [ "16", "Ssid_016", "Ssid_016", "MacN", "1", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "12:15", "37", "0.389", "0.767", "19", "0", "1", "22.50081", "8.048649" ], [ "17", "Ssid_017", "Ssid_017", "MacN", "2", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "12:22", "38", "0.402", "0.865", "25", "0", "1", "22.543476", "8.055289" ], [ "18", "Ssid_018", "Ssid_018", "MacN", "4", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "12:23", "40", "0.299", "0.726", "18", "1", "0", "22.54923", "8.053229" ], [ "19", "Ssid_019", "Ssid_019", "MacN", "5", "1", "3", "controlpH", "controltemp", "controlpH_controltemp", "CC", "12:22", "41", "0.445", "0.879", "25", "0", "1", "22.543476", "8.055289" ], [ "20", "Ssid_020", "Ssid_020", "Emerald", "1", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:31", "22", "0.279", "0.607", "21", "0", "1", "32.973378", "8.049314" ], [ "21", "Ssid_021", "Ssid_021", "Emerald", "2", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:32", "23", "0.206", "0.414", "16", "0", "0.4", "32.988851", "8.04916" ], [ "22", "Ssid_022", "Ssid_022", "Emerald", "3", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:34", "24", "0.301", "0.518", "16", "0", "0", "32.98118", "8.049152" ], [ "23", "Ssid_023", "Ssid_023", "Emerald", "4", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:34", "25", "0.195", "0.436", "21", "0", "1", "32.98118", "8.049152" ], [ "24", "Ssid_024", "Ssid_024", "Emerald", "5", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:36", "26", "0.281", "0.508", "21", "0", "1", "32.996735", "8.050692" ], [ "25", "Ssid_025", "Ssid_025", "Rainbow", "1", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:38", "27", "0.334", "0.641", "13", "0", "1", "32.997342", "8.051248" ], [ "26", "Ssid_026", "Ssid_026", "Rainbow", "2", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:39", "28", "0.359", "0.677", "9", "1", "0", "33.023142", "8.050659" ], [ "27", "Ssid_027", "Ssid_027", "Rainbow", "3", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:41", "29", "0.284", "0.973", "7", "0", "0", "33.013946", "8.05067" ], [ "28", "Ssid_028", "Ssid_028", "Rainbow", "4", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:41", "30", "0.269", "0.51", "21", "0", "1", "33.013946", "8.05067" ], [ "29", "Ssid_029", "Ssid_029", "Rainbow", "5", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:43", "31", "0.262", "0.503", "21", "0", "1", "32.997735", "8.050595" ], [ "30", "Ssid_030", "Ssid_030", "Star", "1", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:47", "32", "0.278", "0.563", "21", "0", "1", "32.964428", "8.05203" ], [ "31", "Ssid_031", "Ssid_031", "Star", "2", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:46", "33", "0.297", "0.569", "21", "0", "1", "32.979934", "8.052668" ], [ "32", "Ssid_032", "Ssid_032", "Star", "3", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:49", "34", "0.308", "0.632", "21", "0", "1", "32.949266", "8.05349" ], [ "33", "Ssid_033", "Ssid_033", "Star", "4", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:49", "35", "0.217", "0.446", "21", "0", "1", "32.949266", "8.05349" ], [ "34", "Ssid_034", "Ssid_034", "Star", "5", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:51", "36", "0.369", "0.665", "30", "0", "1", "32.960167", "8.052777" ], [ "35", "Ssid_035", "Ssid_035", "MacN", "1", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:51", "37", "0.255", "0.47", "11", "0", "1", "32.960167", "8.052777" ], [ "36", "Ssid_036", "Ssid_036", "MacN", "2", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:53", "38", "0.283", "0.56", "16", "0", "1", "32.993998", "8.052947" ], [ "37", "Ssid_037", "Ssid_037", "MacN", "3", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:55", "39", "0.292", "0.615", "21", "0", "0", "32.989228", "8.053293" ], [ "38", "Ssid_038", "Ssid_038", "MacN", "4", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:56", "40", "0.233", "0.502", "29", "0.5", "0.5", "32.996998", "8.051359" ], [ "39", "Ssid_039", "Ssid_039", "MacN", "5", "2", "5", "controlpH", "hightemp", "controlpH_hightemp", "CH", "12:57", "41", "0.338", "0.634", "28", "0", "1", "33.004685", "8.049731" ], [ "40", "Ssid_040", "Ssid_040", "Emerald", "1", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:02", "22", "0.316", "0.635", "25", "0", "1", "22.509874", "7.803966" ], [ "41", "Ssid_041", "Ssid_041", "Emerald", "2", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:03", "23", "0.296", "0.605", "12", "0", "0", "22.521856", "7.803734" ], [ "42", "Ssid_042", "Ssid_042", "Emerald", "3", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:04", "24", "0.373", "0.693", "12", "0", "0.9", "22.53256", "7.803478" ], [ "43", "Ssid_043", "Ssid_043", "Emerald", "4", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:06", "25", "0.345", "0.639", "25", "0", "0.98", "22.556229", "7.806334" ], [ "44", "Ssid_044", "Ssid_044", "Emerald", "5", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:07", "26", "0.376", "0.667", "12", "0", "0", "22.556212", "7.805907" ], [ "45", "Ssid_045", "Ssid_045", "Rainbow", "1", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:10", "27", "0.31", "0.613", "17", "0", "1", "22.583783", "7.809321" ], [ "46", "Ssid_046", "Ssid_046", "Rainbow", "2", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:11", "28", "0.333", "0.625", "12", "1", "0", "22.59047", "7.81134" ], [ "47", "Ssid_047", "Ssid_047", "Rainbow", "3", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:12", "29", "0.268", "0.751", "3", "0", "0", "22.606288", "7.813672" ], [ "48", "Ssid_048", "Ssid_048", "Rainbow", "4", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:14", "30", "0.312", "0.591", "25", "0", "1", "22.63212", "7.816172" ], [ "49", "Ssid_049", "Ssid_049", "Rainbow", "5", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:16", "31", "0.294", "0.536", "25", "0", "1", "22.598649", "7.816183" ], [ "50", "Ssid_050", "Ssid_050", "Star", "1", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:21", "32", "0.345", "0.682", "25", "0", "1", "22.456144", "7.813335" ], [ "51", "Ssid_051", "Ssid_051", "Star", "2", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:23", "33", "0.367", "0.707", "11", "0", "1", "22.435097", "7.811994" ], [ "52", "Ssid_052", "Ssid_052", "Star", "3", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:24", "34", "0.324", "0.64", "25", "0", "1", "22.437146", "7.808319" ], [ "53", "Ssid_053", "Ssid_053", "Star", "4", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:25", "35", "0.345", "0.69", "11", "0", "1", "22.430131", "7.807912" ], [ "54", "Ssid_054", "Ssid_054", "Star", "5", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:26", "36", "0.433", "0.752", "25", "0", "1", "22.441637", "7.808838" ], [ "55", "Ssid_055", "Ssid_055", "MacN", "1", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:27", "37", "0.341", "0.628", "17", "0", "1", "22.447702", "7.810189" ], [ "56", "Ssid_056", "Ssid_056", "MacN", "2", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:29", "38", "0.359", "0.731", "25", "0", "1", "22.459897", "7.811078" ], [ "57", "Ssid_057", "Ssid_057", "MacN", "3", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:30", "39", "0.348", "0.678", "12", "0", "1", "22.460389", "7.810483" ], [ "58", "Ssid_058", "Ssid_058", "MacN", "4", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:30", "40", "0.231", "0.519", "27", "0.88", "0.12", "22.460389", "7.810483" ], [ "59", "Ssid_059", "Ssid_059", "MacN", "5", "3", "7", "lowpH", "controltemp", "lowpH_controltemp", "LC", "1:30", "41", "0.38", "0.735", "25", "0", "1", "22.460389", "7.810483" ], [ "60", "Ssid_060", "Ssid_060", "Emerald", "1", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:35", "22", "0.283", "0.597", "11", "0", "1", "32.941185", "7.809012" ], [ "61", "Ssid_061", "Ssid_061", "Emerald", "2", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:36", "23", "0.197", "0.421", "16", "0", "0", "32.954266", "7.809214" ], [ "62", "Ssid_062", "Ssid_062", "Emerald", "3", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:37", "24", "0.274", "0.486", "11", "0", "0", "32.949545", "7.809134" ], [ "63", "Ssid_063", "Ssid_063", "Emerald", "4", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:38", "25", "0.289", "0.568", "21", "0", "1", "32.962183", "7.810044" ], [ "64", "Ssid_064", "Ssid_064", "Emerald", "5", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:39", "26", "0.299", "0.568", "11", "0", "0.81", "32.968805", "7.810577" ], [ "65", "Ssid_065", "Ssid_065", "Rainbow", "1", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:42", "27", "0.322", "0.635", "11", "0", "1", "32.985425", "7.810489" ], [ "66", "Ssid_066", "Ssid_066", "Rainbow", "2", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:43", "28", "0.279", "0.521", "7", "1", "0", "33.01134", "7.811204" ], [ "67", "Ssid_067", "Ssid_067", "Rainbow", "3", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:44", "29", "0.271", "0.79", "7", "0", "1", "33.015061", "7.810553" ], [ "68", "Ssid_068", "Ssid_068", "Rainbow", "4", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:45", "30", "0.284", "0.545", "21", "0", "1", "33.02047", "7.810553" ], [ "69", "Ssid_069", "Ssid_069", "Rainbow", "5", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:46", "31", "0.281", "0.52", "21", "0", "1", "33.019503", "7.811623" ], [ "70", "Ssid_070", "Ssid_070", "Star", "1", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:50", "32", "0.272", "0.529", "21", "0", "0", "32.97859", "7.810657" ], [ "71", "Ssid_071", "Ssid_071", "Star", "2", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:51", "33", "0.291", "0.559", "21", "0", "1", "32.96987", "7.809632" ], [ "72", "Ssid_072", "Ssid_072", "Star", "3", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:52", "34", "0.254", "0.525", "21", "0", "1", "32.947414", "7.786134" ], [ "73", "Ssid_073", "Ssid_073", "Star", "4", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:54", "35", "0.222", "0.455", "23", "0", "1", "32.956462", "7.809606" ], [ "74", "Ssid_074", "Ssid_074", "Star", "5", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:55", "36", "0.35", "0.627", "30", "0", "1", "32.953462", "7.806517" ], [ "75", "Ssid_075", "Ssid_075", "MacN", "1", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:57", "37", "0.261", "0.492", "11", "0", "1", "32.965281", "7.712204" ], [ "76", "Ssid_076", "Ssid_076", "MacN", "2", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "1:59", "38", "0.219", "0.45", "16", "0", "0", "32.997162", "7.742619" ], [ "77", "Ssid_077", "Ssid_077", "MacN", "3", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "2:00", "39", "0.301", "0.604", "21", "0", "1", "32.997391", "7.771193" ], [ "78", "Ssid_078", "Ssid_078", "MacN", "4", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "2:01", "40", "0.245", "0.562", "30", "1", "0", "33.014586", "7.812524" ], [ "79", "Ssid_079", "Ssid_079", "MacN", "5", "4", "2", "lowpH", "hightemp", "lowpH_hightemp", "LH", "2:02", "41", "0.331", "0.649", "23", "0", "0", "33.00339", "7.814244" ] ];
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
