# >> Experiment files
scenario = "Tinnits_calibration";
pcl_file = "Tinnitus_calibration.pcl";
# << Experiment files

# >> Experiment; config
scenario_type = trials;     
response_logging = log_active;  
response_matching = simple_matching; 
randomize_trials = true;
active_buttons = 6;
button_codes =253, 254, 255, 253, 254, 255;
default_output_port = 1; 
response_port_output = true; 
write_codes = true; 
pulse_width = 10; 
default_background_color = 0,0,0;   
default_font_size = 15;
default_text_color = 255,255,255;
default_font = "Tahoma";      
# << Experiment; config     

                            
# >>Init Experiment *************************************
begin;

#picture {} default;
picture {text {caption = " "; font_size=40; font_color = 255,255,255;}; x=0;y=0;} default;
picture {text {caption = "x"; font_size=40; font_color = 255,255,255;}; x=0;y=0;} afixX;
picture {text {caption = "The test is now finished. Thank you"; font_size=40; font_color = 255,255,255;}; x=0;y=0;} endOfTest;
picture {text {caption = "How many deviating sounds were there?"; font_size=40; font_color = 255,255,255;}; x=0;y=0;} countSounds;
picture {text {caption = "My test"; font_size=40; font_color = 255,255,255;}; x=0;y=0;} mytest;
picture {text {caption = "My end"; font_size=40; font_color = 255,255,255;}; x=0;y=0;} mytest2;

text { caption = "Hello!"; preload = false; font_size=24;} myStatusText;

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

trial {
   trial_duration =1;
   picture {
      text myStatusText;
		x = 0; y = 0;
   } statusPic;
} tr_statusText;

trial {
	all_responses = true;
	trial_duration = forever;
	trial_type = first_response; 
	picture {
		text {caption = "The test is now finished. Thank you"; font_size=42;}; x = 0; y = 0;
	};         
	time = 0;
} tr_endoftest;




# OLD STUFF

#ITI
trial {
	all_responses = false;
	trial_duration = 1;
	trial_type = first_response; 
	picture {text {caption = " "; font_size=40; font_color = 255,255,255;}; x=0;y=0;} blank;
} tr_iti;

#FIXATION
trial {
	all_responses = false;
	trial_duration = 1;
	trial_type = first_response; 
	stimulus_event{ 
		picture {text {caption = "LoadSound"; font_size=40; font_color = 255,100,100;}; x=0;y=0;
		} afix1; 
		code = "LoadSound";
		port_code = 64;
	} ev_loadSound; 
} tr_loadSound;

#Main trial
trial {     
	all_responses = false;
	trial_type = first_response;
	trial_duration = 50;	
   stimulus_event{ 
      picture {
         text {caption = "Standard sound"; font_size=40; font_color = 100,255,100;}; x = 0; y = 0;
         } pic_Stimulation;  
      code = "Standard sound";
   } ev_stimulation; 
} tr_stimulation; 

trial {     
	all_responses = false;
	trial_type = first_response;
	trial_duration = 1;	
   stimulus_event{ 
      picture {
         text {caption = "Oddball"; font_size=40; font_color = 100,255,100;}; x = 0; y = 0;
         } pic_Oddball;  
      code = "Oddball";
      port_code = 2;
   } ev_oddball; 
} tr_oddball; 

#Feedback break
trial {
	all_responses = true;
	trial_duration = forever;
	trial_type = first_response; 
	picture {
		text {caption = "Check how tired / alert the subject is"; font_size=42;}; x = 0; y = 0;
	};         
	time = 0;
} tr_break;