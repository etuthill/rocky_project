// Motor test and print calibration data to serial monitor
// To run this demo, you will need to install the LSM6 library:
// https://github.com/pololu/lsm6-arduino

#include <Balboa32U4.h>
#include <Wire.h>
#include <LSM6.h>

LSM6 imu;
Balboa32U4Motors motors;
Balboa32U4Encoders encoders;

const uint8_t UPDATE_TIME_MS = 10;
const uint16_t SQUARE_WAVE_PERIOD_MS = 2000;
const uint8_t NUM_PERIODS = 20;

//time variables
uint32_t startTime_ms;
uint32_t prevTime_ms;
uint32_t curTime_ms;
float curTime_s;

//encoder count variables
int32_t countsLeft = 0;
int32_t countsRight = 0;
int32_t deltaCountsLeft;
int32_t deltaCountsRight;

//distance travelled variables
float distLeft_m;
float distRight_m;
float deltaDistLeft_m;
float deltaDistRight_m;

//speed variables
float speedLeft_mps;
float speedRight_mps;

//speed command variables
int16_t leftSpeedCommand;
int16_t rightSpeedCommand;


void setup()
{
  Serial.begin(9600);
  startTime_ms = millis();
  prevTime_ms = startTime_ms;
  delay(2+SQUARE_WAVE_PERIOD_MS/2);
}

void loop()
{

  curTime_ms = millis();                   // get the current time in miliseconds
  curTime_s = ((float) curTime_ms) / ((float) 1000.0);

  //alternate motors between turning on and off
  if (curTime_ms - startTime_ms > NUM_PERIODS*SQUARE_WAVE_PERIOD_MS){
    leftSpeedCommand = 0;
    rightSpeedCommand = 0;
  }
  else{
    if (2*((curTime_ms - startTime_ms)%SQUARE_WAVE_PERIOD_MS)<=SQUARE_WAVE_PERIOD_MS){
      leftSpeedCommand = 300;
      rightSpeedCommand = 300;
    }
    else{
      leftSpeedCommand = 0;
      rightSpeedCommand = 0; 
    }
    
  }
  motors.setSpeeds(leftSpeedCommand, rightSpeedCommand);

  if (curTime_ms - prevTime_ms >= UPDATE_TIME_MS  && curTime_ms - startTime_ms < NUM_PERIODS*SQUARE_WAVE_PERIOD_MS) {
    integrateEncoders();
    motorUnitConversion();
    prevTime_ms = curTime_ms;

    Serial.print(curTime_s,3); //print time (in seconds)
    Serial.print(" ");
    Serial.print(leftSpeedCommand); //print left wheel commanded speed (dimensionless)
    Serial.print(" ");
    Serial.print(speedLeft_mps,6);  //print left wheel speed (in meters/second)
    Serial.print(" ");
    Serial.print(rightSpeedCommand); //print right wheel commanded speed (dimensionless)
    Serial.print(" ");
    Serial.print(speedRight_mps,6); //print right wheel speed (in meters/second)
    Serial.println("");

  }
}

void motorUnitConversion()
{
  //wheel radius = 8 cm = .08 meters
  //motor gearbox ratio = 38437/507 approx 75/1
  //external gearbox ratio = 45/21
  //combined gearbox ratio approx 162.5/1
  //number of encoder counts per motor rotation = 12
  //distance/#ticks = pi*radius / (combined ratio* ticks/rotation)
  //distance/#ticks = 1.28921060e-4
  
  //total distance travelled by each wheel
  distLeft_m = ((float)countsLeft) * ((float) 1.289210607555819) / 10000.0;
  distRight_m = ((float)countsRight) * ((float) 1.289210607555819) / 10000.0;

  //incremental distance travelled by each wheel since last timestep
  deltaDistLeft_m = ((float)deltaCountsLeft) * ((float) 1.289210607555819) / 10000.0;
  deltaDistRight_m = ((float)deltaCountsRight) * ((float) 1.289210607555819) / 10000.0;

  //compute speed in meters per second (note that time is in ms)
  speedLeft_mps = deltaDistLeft_m * ((float) 1000.0) /((float) (curTime_ms - prevTime_ms));
  speedRight_mps = deltaDistRight_m * ((float) 1000.0) /((float) (curTime_ms - prevTime_ms));
}

void integrateEncoders()
{
  deltaCountsLeft = encoders.getCountsAndResetLeft();
  deltaCountsRight = encoders.getCountsAndResetRight();

  countsLeft += deltaCountsLeft;
  countsRight += deltaCountsRight; 
}