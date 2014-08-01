/*	This code is based on a number of papers on battery lifetime modelling.
	Each function that implements an equation taken from those papers is
	referenced for easier understanding. */

#include "Timer.h"
#include <SPI.h>
#include <SD.h>
#include <stdlib.h>
#include <math.h>

/********* DATA OBTAINED FROM DATASHEET *********/
float nominal_cell_voltage = 2.23;
int nominal_battery_temperature = 298; // Kelvin
float nominal_capacity = 7; // Ah
int cell_number = 6;
float battery_height = 0.3; // in meters
float iec_cycles = 420;
float float_lifetime = 5*24*365; // in hours

//*********** DATA NOT OBTAINED FROM DATASHEET ***********/

float corr_voltage_0 = 1.75; // For 30% acid concentration
float internal_resistance_c = 0.42; // [Ohms Ah]
float internal_resistance_d = 0.699; // [Ohms Ah]
float charge_ov_coefficient = 0.888; // Non-dimensional
float discharge_ov_coefficient = 0.0464; // Non-dimensional
float normalized_capacity_c = 1.001; // Non-dimensional
float normalized_capacity_d = 1.75; // Non-dimensional
float oc_voltage_gradient = 0.076; // Volts


/*********** CHANGEABLE PARAMETERS ***********/
float c_deg_lim = 0.8;
float corrosion_time_step = 0.016667; // [hours]
float degradation_time_step = 0.016667; // [hours]

/******* PARAMETERS FOR I_gas CALCULATION *******/

// Normalized Gassing Current for a 100Ah VRLA battery in Amps
float i_gas_0 = 0.01;
int voltage_coefficient = 11; // 1/V
float temperature_coefficient = 0.06; // 1/Kelvin


/*******  GLOBAL VARIABLES USED IN CODE *******/
float soc_limit = 0.9;

// Measurement Variables
float battery_voltage, battery_temperature, battery_current, ambient_temperature, soc;
float current_time_step = 0.25; // In Seconds
float *current_vector; float current_vector_size = 240; int current_vector_pos = 0;
int soc_status[2] = {0,0}; // Helps us determine changes in battery charge state

// C_deg Variables
float c_deg, current_factor, f_soc, f_plus, f_minus, f_acid, soc_min, soc_max,
		f_strat, z_w, time_since_last_charge, n_bad_charges;


// C_corr Variables
float corrosion_temperature_0 = 298; // Corrosion reference temperature in K
float dw, dw_lim; // Initial value for corrosion layer growth
float c_corr;

// C_remaining Variables
float c_remaining;
float c_remaining_history[2] = {1.75,1.75};


// Arduino Variables
const int chipSelect = 53;
Timer t;
File dataFile;


// TODO: Review function calls and eliminate unecessary parameters
float integrate(float *y, float time_step, int sample_size) {
	// Integrates using 3/8 Simpson method
	float sum=0,h,temp;
	int i,j,k=0,l=0;
	
	sum = sum + y[0];

	for(i=1;i<(sample_size-1);i++) {
		if(k==0 || l==0) {
			sum = sum + 3 * y[i];
			if(k==1) l=1;
			k=1;
		}
		else {
			sum = sum + 2 * y[i];
			k=0;
			l=0;
		}
	}
	sum = sum + y[i];
	sum = sum * (3*time_step/8);
	return sum;
}

float measure_battery_current() {
	float aggregate_data[100], data_sum = 0, current_constant = 0.064, current_offset = 2.5075;
	int i;
	for (i=0;i<100;i++) {
		aggregate_data[i] = analogRead(A3);
		data_sum += aggregate_data[i];
	}

	float current_value = data_sum/100;
	current_value = (current_offset - current_value * (5.00 / 1023.00))/current_constant;
	if (fabs(current_value) < 0.1) current_value = 0;

	log_data();
	return current_value;
}

float measure_battery_voltage() {
	float aggregate_data[100], data_sum = 0, voltage_constant = 2.97;
	int i;
	for (i=0;i<100;i++) {
		aggregate_data[i] = analogRead(A0);
		data_sum += aggregate_data[i];
	}
	float voltage_value = data_sum/100;
	return (voltage_value * (5.00 / 1023.00)) * voltage_constant;
}

float measure_ambient_temperature() {
	float reading = analogRead(A1);
	reading = (5.0 * reading * 100.0)/1024.0;
	return reading;
}

float measure_battery_temperature() {
	float reading = analogRead(A2);
	reading = (5.0 * reading * 100.0)/1024.0;
	return reading;
}

/* Schiffer et al. (2007) / Section 3.2 / Equation #6 */
float calculate_i_gas(float battery_voltage, float battery_temperature) {
	float cell_voltage = battery_voltage / cell_number;
	
	float gassing_current = (nominal_capacity/100)*i_gas_0*
		exp(voltage_coefficient*(cell_voltage - nominal_cell_voltage) + 
			temperature_coefficient*(battery_temperature - nominal_battery_temperature));

	return gassing_current;
}

/* Dufo-López et al. (2014) / Section 4.3.1.1 / Equation #13 */
void calculate_n_bad_charges() {
	n_bad_charges += (0.0025 - pow((0.95 - soc_max),2))/0.0025;
	if (soc > 0.9999) n_bad_charges = 0;
}

/* Dufo-López et al. (2014) / Section 4.3.1.1 / Equation #12 */
void calculate_current_factor() {
	calculate_n_bad_charges();
	current_factor = sqrt((nominal_capacity/10)/fabs(battery_current)) * 
		exp(n_bad_charges/(3*3.6));
}

/* Schiffer et al. (2007) / Section 3.2 / Equation #7 */
void calculate_soc(float *current_values, float current_time_step, int sample_size) {
	// Positive currents will charge the battery
	float i_gas = calculate_i_gas(battery_voltage, battery_temperature);
	int i;

	for (i=0;i<sample_size;i++) {
		current_values[i] -= i_gas;
	}

	soc += integrate(current_values, current_time_step, sample_size)/(nominal_capacity*3600); // *3600 to convert Ah to A*seconds
	
	/* TODO: MAYBE REMOVE THIS, IMPLEMENT CLOCK */
	time_since_last_charge += 1/60;
	
	// Cap SOC at 1
	if (soc > 1) soc = 1;

	if (soc >= soc_limit) {
		time_since_last_charge = 0;

		// Battery is now charged, must update the status
		soc_status[0] = soc_status[1];
		soc_status[1] = 1;

		if (soc_status[0] == 0 && soc_status[1] == 1) soc_max = 0; // resetting soc_max
		if (soc > soc_max) soc_max = soc;
	}
	else {
		
		// Battery is not charged, must update the status
		soc_status[0] = soc_status[1];
		soc_status[1] = 0;
		
		if (soc_status[0] == 1 && soc_status[1] == 0) { // Detects change in charged state
			soc_min = 1;
			calculate_current_factor();
			if (soc < soc_min) soc_min = soc;
		}
	}
}

/* Dufo-López et al. (2014) / Section 4.3.1.1 / Equation #11 */
void calculate_f_soc() {
	float c_soc_0 = 0.00006614; // hour^-1
	float c_soc_min = 0.003307; // hour^-1

	f_soc = (1 + (c_soc_0 + c_soc_min * (1 - soc_min) * current_factor * time_since_last_charge));
}

/* Schiffer et al. (2007) / Section 3.4.3.1 / Equation #23 */
void calculate_f_plus() {
	// Will return 0 on the charging process, do the math otherwise
	float c_plus = (1.0/30.0);
	if (battery_current >= 0) f_plus = 0;
	else {
		f_plus = (c_plus * (1 - soc_min) * exp((-3) * f_strat) * (fabs(battery_current)/(nominal_capacity/10)));
	}
}

/* Schiffer et al. (2007) / Section 3.4.3.2 / Equations #24, 27 & 28 */
void calculate_f_minus() {
	float c_minus = 0.1 , voltage_constant = 0.167;
	double stratification_removal = pow(10,-9);
	float cell_voltage = battery_voltage/cell_number;
	float battery_temperature = battery_temperature;
	
	//TODO: Must investigate the multiplication by 1 to take into account battery ageing
	float f_minus_gassing = c_minus * sqrt(100/nominal_capacity) * 1 * 
		exp(voltage_constant * (cell_voltage - 2.5) + temperature_coefficient * 
		(( battery_temperature + 273) - nominal_battery_temperature));
	double f_minus_diffusion = 8.0 * (stratification_removal / pow(battery_height,2)) * 
		f_strat * pow(2,((battery_temperature - 20)/10.0));
	f_minus = (f_minus_gassing + (float)f_minus_diffusion); // TODO: Fix this crappy casting
}

/* Schiffer et al. (2007) / Section 3.4.3 / Equation #21 */
void calculate_f_strat(float time_step) { // time_step must be given in hours
	calculate_f_plus(); calculate_f_minus();
	f_strat += ( f_plus - f_minus) * time_step;
	if (f_strat < 0) f_strat = 0;
}

/* Schiffer et al. (2007) / Section 3.4.3 / Equation #22 */
void calculate_f_acid(float time_step) { // time_step must be given in hours
	calculate_f_strat(time_step);
	f_acid = 1 +  f_strat * sqrt( (nominal_capacity/10)/fabs(battery_current) );
	if (battery_current == 0) f_acid = 1; // TODO: This is unreviewed policy, must be validated
}

/* Dufo-López et al. (2014) / Section 4.3.1 / Equation #10 */
void calculate_z_w(float time_step) { // time_step must be given in hours
	calculate_f_soc();
	calculate_f_acid(time_step);
	z_w += (fabs(battery_current) * f_soc * f_acid * time_step) / nominal_capacity;
}

/* Dufo-López et al. (2014) / Section 4.3.1 / Equation #9 */
void calculate_c_deg() { // time_step must be given in hours
	calculate_z_w(degradation_time_step);
	c_deg = (c_deg_lim * exp((-5)*(1 - (z_w/(1.6 * iec_cycles)))));
}

/* Schiffer et al. (2007) / Section 3.1 / Equation #5 */
float calculate_corr_voltage() {
	/* TODO: Investigate this function's return values. For high SOCs, it
	   returns crazy high corrosion voltages */
	float corr_voltage = 0;
	if (battery_current > 0) {
		corr_voltage = corr_voltage_0 - (10.0/13.0) * oc_voltage_gradient * (1 - soc) +
		 0.5 * internal_resistance_c * (battery_current / nominal_capacity) +
		 0.5 * internal_resistance_c * charge_ov_coefficient *
		 (battery_current / nominal_capacity) * (soc / (normalized_capacity_c - soc));
	} 
	else {
		corr_voltage = corr_voltage_0 - (10.0/13.0) * oc_voltage_gradient * (1 - soc) +
		 0.5 * internal_resistance_d * (battery_current / nominal_capacity) +
		 0.5 * internal_resistance_d * discharge_ov_coefficient *
		 (battery_current / nominal_capacity) * ((1 - soc) / (normalized_capacity_d - (1 - soc)));
	}
	return corr_voltage;
}

/* This is a manual fit of the Lander Curve. */
/* Bindner et al. (2005) / Section 5.2.1.1 / Figure #11 */
float calculate_corrosion_speed(float corrosion_voltage) {
	if (corrosion_voltage < 0) return 0;

	float lander_curve_x[13] = {0, 0.5, 1.576, 1.598, 1.664, 1.686, 1.709, 1.739, 1.950, 
		2.029, 2.099, 2.151, 2.184};
	float lander_curve_y[13] = {0.053305, 0.055084, 0.346755, 0.666663, 0.368389, 0.347145, 
		0.048716, 0.016842, 0.017591, 0.177787, 0.615137, 2.00125, 5.00777};

	
	int i = 0;
	while(corrosion_voltage >= lander_curve_x[i]) i++;
	if (i>=13) return 5;
	float slope = (lander_curve_y[i] - lander_curve_y[i-1]) / (lander_curve_x[i] - lander_curve_x[i-1]);

	return (lander_curve_y[i-1] + slope*(corrosion_voltage - lander_curve_x[i-1]));
}

/* Schiffer et al. (2007) / Section 3.3 / Equation #9 */
float calculate_k_s(float corrosion_voltage) {
	float k_s_temperature = log(2)/15; // Corrosion speed doubles with a 15K temperature increase
	float k_s = calculate_corrosion_speed(corrosion_voltage) * 
		exp(k_s_temperature * (battery_temperature + 273 - corrosion_temperature_0));
	
	return k_s;
}

/* Schiffer et al. (2007) / Section 3.3 / Equation #8 */
void calculate_dw() {
	float corrosion_voltage = calculate_corr_voltage();
	float k_s = calculate_k_s (corrosion_voltage);
	
	if (corrosion_voltage < 1.74)
		dw = k_s * pow((pow(dw/k_s,(1.0/0.6)) + corrosion_time_step), 0.6);
	else 
		dw = dw + k_s * corrosion_time_step;
}

/* Dufo-López et al. (2014) / Section 4.3.2 / Equation #16 */
void calculate_c_corr() {
	calculate_dw();
	c_corr = (0.2 * dw/dw_lim);
}

void calculate_remaining_lifetime() {
	/* TODO: Improve this remaining lifetime algorythm. So far, it assumes
	   a linear decrease in c_remaining. Maybe project usage patterns? */

	float slope = (c_remaining_history[1] - c_remaining_history[0])/ corrosion_time_step;

	float remaining_lifetime = (0.8 - 1.75) / (slope * 24 *365);
	printf("%.2f\n",remaining_lifetime);
}

/* Schiffer et al. (2007) / Section 3.5 / Equation #30 */
void calculate_c_remaining(void *context) {
	calculate_c_corr(); calculate_c_deg();
	c_remaining =  1.75 -  c_corr - c_deg;

	// This is used for remaining lifetime prediction
	c_remaining_history[0] = c_remaining_history[1];
	c_remaining_history[1] = c_remaining;
	calculate_remaining_lifetime();
}

void take_measurements(void *context) {
	battery_voltage = measure_battery_voltage();
	battery_temperature = measure_battery_temperature();
	ambient_temperature = measure_ambient_temperature();
	battery_current = measure_battery_current();

	if (current_vector_pos == current_vector_size) {
		calculate_soc(current_vector, current_time_step, current_vector_size);
		current_vector_pos = 0;
	}
	else {
		current_vector[current_vector_pos] = battery_current;
		current_vector_pos++;
	}
}
		
float read_value() {
	char buffer[10]; 
	int buffer_pos = 0;
	char read_data;
	float final_value;

	read_data = dataFile.read();

	while ( (read_data != ',') && (read_data != '\n') ) {
		if (read_data == ' ');
		else {
			buffer[buffer_pos] = read_data;
			buffer_pos++;
		}
		read_data = dataFile.read();
	}
	
	final_value = atof(buffer);
	return final_value;
}

int rewind_line() {
	int number_columns = 0;

	long dataPos = dataFile.position();
	dataPos -= 2;
	dataFile.seek(dataPos);
	char read_data = dataFile.read();

	while(read_data != '\n') {
		dataPos -= 1;
		dataFile.seek(dataPos);
		if (read_data == ',') number_columns++;

		read_data = dataFile.read();
	}
	return (number_columns + 1);
}

/* TODO: DEVELOP INITIAL SOC MEASUREMENT */
/* TODO: GATHER DATA BEFORE RESET FROM TEXT FILE */
void initial_setup() {
	delay(5000);

	n_bad_charges = 0; soc_min = 1; soc_max = 0;
	f_strat = 0; z_w = 0; time_since_last_charge = 0;

	current_vector = (float *)malloc(current_vector_size * sizeof(float));

	/* Bindner et al. (2005) / Section 5.3.7.2 / Equation #32 */
	dw_lim = float_lifetime * calculate_k_s(nominal_cell_voltage);

	// Equation obtained from linear regression from data relating Voltage and SOC
	battery_voltage = measure_battery_voltage();
	soc = 0.6839*battery_voltage - 7.9077;
	if (soc > 1) soc = 1;
	
	// Gathering previous data from SD card
	if (dataFile.position() > 5) { // Verifies if data file actually has any data, 5 is arbitrary
		Serial.println("Reading Previous Values");	
		int column_read_count = 0;
		float read_values[100];

		int number_columns = rewind_line();
		Serial.println(number_columns);

		while (column_read_count < number_columns) {
			read_values[column_read_count] = read_value();
			Serial.println(read_values[column_read_count]);
			column_read_count++;
		}
		
		// Set variables to correspondent read values
		battery_voltage = read_values[0];
		battery_current = read_values[1];
		soc = read_values[2];
		battery_temperature = read_values[3];
		ambient_temperature = read_values[4];
		current_factor = read_values[5];
		f_soc = read_values[6];
		f_plus = read_values[7];
		f_minus = read_values[8];
		f_acid = read_values[9];
		soc_min = read_values[10];
		soc_max = read_values[11];
		f_strat = read_values[12];
		z_w = read_values[13];
		time_since_last_charge = read_values[14];
		n_bad_charges = read_values[15];
		dw = read_values[16];
		c_deg = read_values[17];
		c_corr = read_values[18];
		c_remaining = read_values[19];


		Serial.println("Finish Reading");
	}

	battery_temperature = measure_battery_temperature();
	ambient_temperature = measure_ambient_temperature();
	battery_current = measure_battery_current();
}


void setup() {
	Serial.begin(9600);

	// BEGIN Initializing SD Card
	while (!Serial);
	Serial.print("Initializing SD card...");
	pinMode(SS, OUTPUT);

	if (!SD.begin(chipSelect)) {
		Serial.println("Card failed, or not present");
		while (1) ;
	}
	Serial.println("card initialized.");

	dataFile = SD.open("datalog.txt", FILE_WRITE);
	if (! dataFile) {
		Serial.println("error opening datalog.txt");
		while (1) ;
	}
	// END Initializing SD Card

	initial_setup();

	int measurement_tick = t.every(current_time_step*1000, take_measurements, 0);
	int c_remaining_tick = t.every(corrosion_time_step*1000*3600, calculate_c_remaining, 0);
}

void log_data() {
	Serial.println("Logging Data");
	char buffer[11];
	String dataString = "";
	
	dataString += dtostrf(battery_voltage, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(battery_current, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(soc, 7, 3, buffer);
	dataString += String(',');	
	dataString += dtostrf(battery_temperature, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(ambient_temperature, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(current_factor, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(f_soc, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(f_plus, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(f_minus, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(f_acid, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(soc_min, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(soc_max, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(f_strat, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(z_w, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(time_since_last_charge, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(n_bad_charges, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(dw, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(c_deg, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(c_corr, 7, 3, buffer);
	dataString += String(',');
	dataString += dtostrf(c_remaining, 7, 3, buffer);

	
	dataFile.println(dataString);
	Serial.println(dataString);

	dataFile.flush();

}

void loop() {
	t.update();
}