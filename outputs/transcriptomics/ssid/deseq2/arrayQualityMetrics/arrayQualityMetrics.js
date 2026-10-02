// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = ["false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","true","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","true","false","false","false"];
var arrayMetadata    = [{"array":"1","sampleNames":"Ssid_001","id":"Ssid_001","site":"Emerald","geno":"1","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"11:36","genotype":"22","propB":"0","propD":"0.98","temp_actual":"22.576029","ph_actual":"8.041769","habitat":"reef"},{"array":"2","sampleNames":"Ssid_002","id":"Ssid_002","site":"Emerald","geno":"2","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"11:40","genotype":"23","propB":"0","propD":"0","temp_actual":"22.606222","ph_actual":"8.040733","habitat":"reef"},{"array":"3","sampleNames":"Ssid_003","id":"Ssid_003","site":"Emerald","geno":"3","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"11:43","genotype":"24","propB":"0","propD":"0.87","temp_actual":"22.567883","ph_actual":"8.044771","habitat":"reef"},{"array":"4","sampleNames":"Ssid_004","id":"Ssid_004","site":"Emerald","geno":"4","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"11:47","genotype":"25","propB":"0","propD":"0.99","temp_actual":"22.490319","ph_actual":"8.045541","habitat":"reef"},{"array":"5","sampleNames":"Ssid_005","id":"Ssid_005","site":"Emerald","geno":"5","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"11:50","genotype":"26","propB":"0","propD":"0.5","temp_actual":"22.487894","ph_actual":"8.045853","habitat":"reef"},{"array":"6","sampleNames":"Ssid_006","id":"Ssid_006","site":"Rainbow","geno":"1","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"11:55","genotype":"27","propB":"0","propD":"1","temp_actual":"22.529544","ph_actual":"8.046972","habitat":"reef"},{"array":"7","sampleNames":"Ssid_007","id":"Ssid_007","site":"Rainbow","geno":"2","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"11:58","genotype":"28","propB":"1","propD":"0","temp_actual":"22.551016","ph_actual":"8.047554","habitat":"reef"},{"array":"8","sampleNames":"Ssid_008","id":"Ssid_008","site":"Rainbow","geno":"3","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"12:01","genotype":"29","propB":"0","propD":"0","temp_actual":"22.568932","ph_actual":"8.039805","habitat":"reef"},{"array":"9","sampleNames":"Ssid_009","id":"Ssid_009","site":"Rainbow","geno":"4","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"12:05","genotype":"30","propB":"0","propD":"1","temp_actual":"22.602862","ph_actual":"8.046859","habitat":"reef"},{"array":"10","sampleNames":"Ssid_010","id":"Ssid_010","site":"Rainbow","geno":"5","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"12:08","genotype":"31","propB":"0","propD":"1","temp_actual":"22.588241","ph_actual":"8.037466","habitat":"reef"},{"array":"11","sampleNames":"Ssid_011","id":"Ssid_011","site":"Star","geno":"1","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"12:07","genotype":"32","propB":"0","propD":"1","temp_actual":"22.605321","ph_actual":"8.045583","habitat":"urban"},{"array":"12","sampleNames":"Ssid_012","id":"Ssid_012","site":"Star","geno":"2","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"12:13","genotype":"33","propB":"0","propD":"1","temp_actual":"22.499794","ph_actual":"8.055674","habitat":"urban"},{"array":"13","sampleNames":"Ssid_013","id":"Ssid_013","site":"Star","geno":"3","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"12:17","genotype":"34","propB":"0","propD":"1","temp_actual":"22.503498","ph_actual":"8.053427","habitat":"urban"},{"array":"14","sampleNames":"Ssid_014","id":"Ssid_014","site":"Star","geno":"4","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"12:18","genotype":"35","propB":"0","propD":"1","temp_actual":"22.5113","ph_actual":"8.056571","habitat":"urban"},{"array":"15","sampleNames":"Ssid_015","id":"Ssid_015","site":"Star","geno":"5","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"12:19","genotype":"36","propB":"0","propD":"1","temp_actual":"22.514956","ph_actual":"8.054615","habitat":"urban"},{"array":"16","sampleNames":"Ssid_016","id":"Ssid_016","site":"MacN","geno":"1","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"12:15","genotype":"37","propB":"0","propD":"1","temp_actual":"22.50081","ph_actual":"8.048649","habitat":"urban"},{"array":"17","sampleNames":"Ssid_017","id":"Ssid_017","site":"MacN","geno":"2","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"12:22","genotype":"38","propB":"0","propD":"1","temp_actual":"22.543476","ph_actual":"8.055289","habitat":"urban"},{"array":"18","sampleNames":"Ssid_018","id":"Ssid_018","site":"MacN","geno":"4","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"12:23","genotype":"40","propB":"1","propD":"0","temp_actual":"22.54923","ph_actual":"8.053229","habitat":"urban"},{"array":"19","sampleNames":"Ssid_019","id":"Ssid_019","site":"MacN","geno":"5","replicate":"1","tank":"3","ph":"controlpH","temp":"controltemp","treatment":"controlpH_controltemp","treat":"CC","time":"12:22","genotype":"41","propB":"0","propD":"1","temp_actual":"22.543476","ph_actual":"8.055289","habitat":"urban"},{"array":"20","sampleNames":"Ssid_020","id":"Ssid_020","site":"Emerald","geno":"1","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:31","genotype":"22","propB":"0","propD":"1","temp_actual":"32.973378","ph_actual":"8.049314","habitat":"reef"},{"array":"21","sampleNames":"Ssid_021","id":"Ssid_021","site":"Emerald","geno":"2","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:32","genotype":"23","propB":"0","propD":"0.4","temp_actual":"32.988851","ph_actual":"8.04916","habitat":"reef"},{"array":"22","sampleNames":"Ssid_022","id":"Ssid_022","site":"Emerald","geno":"3","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:34","genotype":"24","propB":"0","propD":"0","temp_actual":"32.98118","ph_actual":"8.049152","habitat":"reef"},{"array":"23","sampleNames":"Ssid_023","id":"Ssid_023","site":"Emerald","geno":"4","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:34","genotype":"25","propB":"0","propD":"1","temp_actual":"32.98118","ph_actual":"8.049152","habitat":"reef"},{"array":"24","sampleNames":"Ssid_024","id":"Ssid_024","site":"Emerald","geno":"5","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:36","genotype":"26","propB":"0","propD":"1","temp_actual":"32.996735","ph_actual":"8.050692","habitat":"reef"},{"array":"25","sampleNames":"Ssid_025","id":"Ssid_025","site":"Rainbow","geno":"1","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:38","genotype":"27","propB":"0","propD":"1","temp_actual":"32.997342","ph_actual":"8.051248","habitat":"reef"},{"array":"26","sampleNames":"Ssid_026","id":"Ssid_026","site":"Rainbow","geno":"2","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:39","genotype":"28","propB":"1","propD":"0","temp_actual":"33.023142","ph_actual":"8.050659","habitat":"reef"},{"array":"27","sampleNames":"Ssid_027","id":"Ssid_027","site":"Rainbow","geno":"3","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:41","genotype":"29","propB":"0","propD":"0","temp_actual":"33.013946","ph_actual":"8.05067","habitat":"reef"},{"array":"28","sampleNames":"Ssid_028","id":"Ssid_028","site":"Rainbow","geno":"4","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:41","genotype":"30","propB":"0","propD":"1","temp_actual":"33.013946","ph_actual":"8.05067","habitat":"reef"},{"array":"29","sampleNames":"Ssid_029","id":"Ssid_029","site":"Rainbow","geno":"5","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:43","genotype":"31","propB":"0","propD":"1","temp_actual":"32.997735","ph_actual":"8.050595","habitat":"reef"},{"array":"30","sampleNames":"Ssid_030","id":"Ssid_030","site":"Star","geno":"1","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:47","genotype":"32","propB":"0","propD":"1","temp_actual":"32.964428","ph_actual":"8.05203","habitat":"urban"},{"array":"31","sampleNames":"Ssid_031","id":"Ssid_031","site":"Star","geno":"2","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:46","genotype":"33","propB":"0","propD":"1","temp_actual":"32.979934","ph_actual":"8.052668","habitat":"urban"},{"array":"32","sampleNames":"Ssid_032","id":"Ssid_032","site":"Star","geno":"3","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:49","genotype":"34","propB":"0","propD":"1","temp_actual":"32.949266","ph_actual":"8.05349","habitat":"urban"},{"array":"33","sampleNames":"Ssid_033","id":"Ssid_033","site":"Star","geno":"4","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:49","genotype":"35","propB":"0","propD":"1","temp_actual":"32.949266","ph_actual":"8.05349","habitat":"urban"},{"array":"34","sampleNames":"Ssid_034","id":"Ssid_034","site":"Star","geno":"5","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:51","genotype":"36","propB":"0","propD":"1","temp_actual":"32.960167","ph_actual":"8.052777","habitat":"urban"},{"array":"35","sampleNames":"Ssid_035","id":"Ssid_035","site":"MacN","geno":"1","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:51","genotype":"37","propB":"0","propD":"1","temp_actual":"32.960167","ph_actual":"8.052777","habitat":"urban"},{"array":"36","sampleNames":"Ssid_036","id":"Ssid_036","site":"MacN","geno":"2","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:53","genotype":"38","propB":"0","propD":"1","temp_actual":"32.993998","ph_actual":"8.052947","habitat":"urban"},{"array":"37","sampleNames":"Ssid_037","id":"Ssid_037","site":"MacN","geno":"3","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:55","genotype":"39","propB":"0","propD":"0","temp_actual":"32.989228","ph_actual":"8.053293","habitat":"urban"},{"array":"38","sampleNames":"Ssid_038","id":"Ssid_038","site":"MacN","geno":"4","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:56","genotype":"40","propB":"0.5","propD":"0.5","temp_actual":"32.996998","ph_actual":"8.051359","habitat":"urban"},{"array":"39","sampleNames":"Ssid_039","id":"Ssid_039","site":"MacN","geno":"5","replicate":"2","tank":"5","ph":"controlpH","temp":"hightemp","treatment":"controlpH_hightemp","treat":"CH","time":"12:57","genotype":"41","propB":"0","propD":"1","temp_actual":"33.004685","ph_actual":"8.049731","habitat":"urban"},{"array":"40","sampleNames":"Ssid_040","id":"Ssid_040","site":"Emerald","geno":"1","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:02","genotype":"22","propB":"0","propD":"1","temp_actual":"22.509874","ph_actual":"7.803966","habitat":"reef"},{"array":"41","sampleNames":"Ssid_041","id":"Ssid_041","site":"Emerald","geno":"2","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:03","genotype":"23","propB":"0","propD":"0","temp_actual":"22.521856","ph_actual":"7.803734","habitat":"reef"},{"array":"42","sampleNames":"Ssid_042","id":"Ssid_042","site":"Emerald","geno":"3","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:04","genotype":"24","propB":"0","propD":"0.9","temp_actual":"22.53256","ph_actual":"7.803478","habitat":"reef"},{"array":"43","sampleNames":"Ssid_043","id":"Ssid_043","site":"Emerald","geno":"4","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:06","genotype":"25","propB":"0","propD":"0.98","temp_actual":"22.556229","ph_actual":"7.806334","habitat":"reef"},{"array":"44","sampleNames":"Ssid_044","id":"Ssid_044","site":"Emerald","geno":"5","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:07","genotype":"26","propB":"0","propD":"0","temp_actual":"22.556212","ph_actual":"7.805907","habitat":"reef"},{"array":"45","sampleNames":"Ssid_045","id":"Ssid_045","site":"Rainbow","geno":"1","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:10","genotype":"27","propB":"0","propD":"1","temp_actual":"22.583783","ph_actual":"7.809321","habitat":"reef"},{"array":"46","sampleNames":"Ssid_046","id":"Ssid_046","site":"Rainbow","geno":"2","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:11","genotype":"28","propB":"1","propD":"0","temp_actual":"22.59047","ph_actual":"7.81134","habitat":"reef"},{"array":"47","sampleNames":"Ssid_047","id":"Ssid_047","site":"Rainbow","geno":"3","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:12","genotype":"29","propB":"0","propD":"0","temp_actual":"22.606288","ph_actual":"7.813672","habitat":"reef"},{"array":"48","sampleNames":"Ssid_048","id":"Ssid_048","site":"Rainbow","geno":"4","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:14","genotype":"30","propB":"0","propD":"1","temp_actual":"22.63212","ph_actual":"7.816172","habitat":"reef"},{"array":"49","sampleNames":"Ssid_049","id":"Ssid_049","site":"Rainbow","geno":"5","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:16","genotype":"31","propB":"0","propD":"1","temp_actual":"22.598649","ph_actual":"7.816183","habitat":"reef"},{"array":"50","sampleNames":"Ssid_050","id":"Ssid_050","site":"Star","geno":"1","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:21","genotype":"32","propB":"0","propD":"1","temp_actual":"22.456144","ph_actual":"7.813335","habitat":"urban"},{"array":"51","sampleNames":"Ssid_051","id":"Ssid_051","site":"Star","geno":"2","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:23","genotype":"33","propB":"0","propD":"1","temp_actual":"22.435097","ph_actual":"7.811994","habitat":"urban"},{"array":"52","sampleNames":"Ssid_052","id":"Ssid_052","site":"Star","geno":"3","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:24","genotype":"34","propB":"0","propD":"1","temp_actual":"22.437146","ph_actual":"7.808319","habitat":"urban"},{"array":"53","sampleNames":"Ssid_053","id":"Ssid_053","site":"Star","geno":"4","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:25","genotype":"35","propB":"0","propD":"1","temp_actual":"22.430131","ph_actual":"7.807912","habitat":"urban"},{"array":"54","sampleNames":"Ssid_054","id":"Ssid_054","site":"Star","geno":"5","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:26","genotype":"36","propB":"0","propD":"1","temp_actual":"22.441637","ph_actual":"7.808838","habitat":"urban"},{"array":"55","sampleNames":"Ssid_055","id":"Ssid_055","site":"MacN","geno":"1","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:27","genotype":"37","propB":"0","propD":"1","temp_actual":"22.447702","ph_actual":"7.810189","habitat":"urban"},{"array":"56","sampleNames":"Ssid_056","id":"Ssid_056","site":"MacN","geno":"2","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:29","genotype":"38","propB":"0","propD":"1","temp_actual":"22.459897","ph_actual":"7.811078","habitat":"urban"},{"array":"57","sampleNames":"Ssid_057","id":"Ssid_057","site":"MacN","geno":"3","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:30","genotype":"39","propB":"0","propD":"1","temp_actual":"22.460389","ph_actual":"7.810483","habitat":"urban"},{"array":"58","sampleNames":"Ssid_058","id":"Ssid_058","site":"MacN","geno":"4","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:30","genotype":"40","propB":"0.88","propD":"0.12","temp_actual":"22.460389","ph_actual":"7.810483","habitat":"urban"},{"array":"59","sampleNames":"Ssid_059","id":"Ssid_059","site":"MacN","geno":"5","replicate":"3","tank":"7","ph":"lowpH","temp":"controltemp","treatment":"lowpH_controltemp","treat":"LC","time":"1:30","genotype":"41","propB":"0","propD":"1","temp_actual":"22.460389","ph_actual":"7.810483","habitat":"urban"},{"array":"60","sampleNames":"Ssid_060","id":"Ssid_060","site":"Emerald","geno":"1","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:35","genotype":"22","propB":"0","propD":"1","temp_actual":"32.941185","ph_actual":"7.809012","habitat":"reef"},{"array":"61","sampleNames":"Ssid_061","id":"Ssid_061","site":"Emerald","geno":"2","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:36","genotype":"23","propB":"0","propD":"0","temp_actual":"32.954266","ph_actual":"7.809214","habitat":"reef"},{"array":"62","sampleNames":"Ssid_062","id":"Ssid_062","site":"Emerald","geno":"3","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:37","genotype":"24","propB":"0","propD":"0","temp_actual":"32.949545","ph_actual":"7.809134","habitat":"reef"},{"array":"63","sampleNames":"Ssid_063","id":"Ssid_063","site":"Emerald","geno":"4","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:38","genotype":"25","propB":"0","propD":"1","temp_actual":"32.962183","ph_actual":"7.810044","habitat":"reef"},{"array":"64","sampleNames":"Ssid_064","id":"Ssid_064","site":"Emerald","geno":"5","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:39","genotype":"26","propB":"0","propD":"0.81","temp_actual":"32.968805","ph_actual":"7.810577","habitat":"reef"},{"array":"65","sampleNames":"Ssid_065","id":"Ssid_065","site":"Rainbow","geno":"1","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:42","genotype":"27","propB":"0","propD":"1","temp_actual":"32.985425","ph_actual":"7.810489","habitat":"reef"},{"array":"66","sampleNames":"Ssid_066","id":"Ssid_066","site":"Rainbow","geno":"2","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:43","genotype":"28","propB":"1","propD":"0","temp_actual":"33.01134","ph_actual":"7.811204","habitat":"reef"},{"array":"67","sampleNames":"Ssid_067","id":"Ssid_067","site":"Rainbow","geno":"3","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:44","genotype":"29","propB":"0","propD":"1","temp_actual":"33.015061","ph_actual":"7.810553","habitat":"reef"},{"array":"68","sampleNames":"Ssid_068","id":"Ssid_068","site":"Rainbow","geno":"4","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:45","genotype":"30","propB":"0","propD":"1","temp_actual":"33.02047","ph_actual":"7.810553","habitat":"reef"},{"array":"69","sampleNames":"Ssid_069","id":"Ssid_069","site":"Rainbow","geno":"5","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:46","genotype":"31","propB":"0","propD":"1","temp_actual":"33.019503","ph_actual":"7.811623","habitat":"reef"},{"array":"70","sampleNames":"Ssid_070","id":"Ssid_070","site":"Star","geno":"1","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:50","genotype":"32","propB":"0","propD":"0","temp_actual":"32.97859","ph_actual":"7.810657","habitat":"urban"},{"array":"71","sampleNames":"Ssid_071","id":"Ssid_071","site":"Star","geno":"2","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:51","genotype":"33","propB":"0","propD":"1","temp_actual":"32.96987","ph_actual":"7.809632","habitat":"urban"},{"array":"72","sampleNames":"Ssid_072","id":"Ssid_072","site":"Star","geno":"3","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:52","genotype":"34","propB":"0","propD":"1","temp_actual":"32.947414","ph_actual":"7.786134","habitat":"urban"},{"array":"73","sampleNames":"Ssid_073","id":"Ssid_073","site":"Star","geno":"4","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:54","genotype":"35","propB":"0","propD":"1","temp_actual":"32.956462","ph_actual":"7.809606","habitat":"urban"},{"array":"74","sampleNames":"Ssid_074","id":"Ssid_074","site":"Star","geno":"5","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:55","genotype":"36","propB":"0","propD":"1","temp_actual":"32.953462","ph_actual":"7.806517","habitat":"urban"},{"array":"75","sampleNames":"Ssid_075","id":"Ssid_075","site":"MacN","geno":"1","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:57","genotype":"37","propB":"0","propD":"1","temp_actual":"32.965281","ph_actual":"7.712204","habitat":"urban"},{"array":"76","sampleNames":"Ssid_076","id":"Ssid_076","site":"MacN","geno":"2","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"1:59","genotype":"38","propB":"0","propD":"0","temp_actual":"32.997162","ph_actual":"7.742619","habitat":"urban"},{"array":"77","sampleNames":"Ssid_077","id":"Ssid_077","site":"MacN","geno":"3","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"2:00","genotype":"39","propB":"0","propD":"1","temp_actual":"32.997391","ph_actual":"7.771193","habitat":"urban"},{"array":"78","sampleNames":"Ssid_078","id":"Ssid_078","site":"MacN","geno":"4","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"2:01","genotype":"40","propB":"1","propD":"0","temp_actual":"33.014586","ph_actual":"7.812524","habitat":"urban"},{"array":"79","sampleNames":"Ssid_079","id":"Ssid_079","site":"MacN","geno":"5","replicate":"4","tank":"2","ph":"lowpH","temp":"hightemp","treatment":"lowpH_hightemp","treat":"LH","time":"2:02","genotype":"41","propB":"0","propD":"0","temp_actual":"33.00339","ph_actual":"7.814244","habitat":"urban"}];
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
