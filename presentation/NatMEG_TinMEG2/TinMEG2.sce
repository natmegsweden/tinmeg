# >> Experiment files
scenario = "TinMEG2";
pcl_file = "TinMEG2.pcl";
# << Experiment files

# >> Experiment; config
scenario_type = trials;     
response_logging = log_active;  
response_matching = simple_matching; 
randomize_trials = true;
active_buttons = 9;
button_codes = 1,2,3,4,5,6,7,8,9;
default_output_port = 1; 
response_port_output = true; 
write_codes = true; 
pulse_width = 10; 
default_background_color = 0,0,0;   
default_font_size = 20;
default_text_color = 255,255,255;
default_font = "Tahoma";      
# << Experiment; config     

                            
# >>Init Experiment *************************************
begin;

#picture {} default;
picture {text {caption = " "; font_size=40; font_color = 255,255,255;}; x=0;y=0;} default;
picture {text {caption = "The test is now finished. Thank you"; font_size=40; font_color = 255,255,255;}; x=0;y=0;} endOfTest;
text { caption = "Hello!"; preload = false; font_size=42;} myStatusText;

#Start trial
trial {
	all_responses = true;
	trial_duration = forever;
	trial_type = first_response; 
	picture {
		text {caption = "The experiment will now begin"; font_size=42;}; x = 0; y = 0;
	};         
	time = 0;
} tr_start;

text {font_size = 14;
	text_align = align_left;
	caption =
	"Hur sömnig har du känt dig under de senaste 5 minuterna?\n
	1  -  Extremt pigg\n
	2  -  Mycket pigg\n
	3  -  Pigg\n
	4  -  Ganska pigg\n
	5  -  Varken pigg eller sömnig\n
	6  -  Lätt sömnig\n
	7  -  Sömnig men ej ansträngande vara vaken\n
	8  -  Sömnig och något ansträngande att vara vaken\n
	9  -  Mycket sömnig, mycket ansträngande att vara vaken\n
	      kämpar mot sömnen";} kss_t;

# Trial: kss
trial {
	trial_duration = forever;
	# trial_duration = 2000; # for testing 
   trial_type = specific_response;    	# trial ends when response
	terminator_button = 1,2,3,4,5,6,7,8,9;
   # terminator_button = 4, 5, 6, 7, 8, 9, 10, 11, 12; # non-MEG
	stimulus_event {
		picture {
			text kss_t; x = 0; y = 0;
		} kss_p;
		time = 0;
		code = "kss";
		target_button = 1,2,3,4,5,6,7,8,9;
	};
} kss;

trial {
   trial_duration =1;
   picture {
      text myStatusText;
		x = 0; y = 0;
   } statusPic;
} tr_statusText;

#end screen
trial {
	all_responses = true;
	trial_duration = forever;
	trial_type = first_response; 
	picture {
		text {caption = "The test is now finished. Thank you"; font_size=42;}; x = 0; y = 0;
	};         
	time = 0;
} tr_endoftest;

#Feedback break
trial {
	all_responses = true;
	trial_duration = forever; #DEBUG
	#trial_duration = 10000; #DEBUG
	trial_type = first_response; 
	terminator_button = 1;
	picture {
		text {caption = "Check how tired/alert the subject is \n Spacebar to continue"; font_size=42;}; x = 0; y = 0;
	};         
	time = 0;
} tr_break;