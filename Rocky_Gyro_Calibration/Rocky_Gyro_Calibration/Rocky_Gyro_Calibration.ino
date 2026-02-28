// Motor test and print calibration data to serial monitor
// To run this demo, you will need to install the LSM6 library:
// https://github.com/pololu/lsm6-arduino

#include <Balboa32U4.h>
#include <Wire.h>
#include <LSM6.h>
// #include "MadgwickAHRS.h"

// Madgwick filter;
LSM6 imu;
Balboa32U4Motors motors;
Balboa32U4Encoders encoders;

const uint8_t UPDATE_TIME_MS = 10;

// Take 100 measurements initially to calibrate the gyro.
const uint8_t CALIBRATION_ITERATIONS = 100;

//time variables
uint32_t startTime_ms;
uint32_t prevTime_ms;
uint32_t curTime_ms;
float curTime_s;

int32_t gYZero;
int32_t angleIntegrated; // millidegrees
float angleGravity;
float angleGravityTemp;

int32_t gY; // degrees/s
int32_t aX; 
int32_t aZ; 

float angle_rad; //radians
float alpha = .8;


void setup()
{
  Serial.begin(9600);

  gyroSetup();

  startTime_ms = millis();
  prevTime_ms = startTime_ms;
}

void loop()
{
  // get the current time in miliseconds
  curTime_ms = millis();        
  curTime_s = ((float) curTime_ms) / ((float) 1000.0);

  if (curTime_ms - prevTime_ms >= UPDATE_TIME_MS) {
    imu.read();
    integrateGyro();

    Serial.print(curTime_s,3); //print time (in seconds)
    Serial.print(" ");
 
    Serial.print(angle_rad,6);
    Serial.println("");

    prevTime_ms = curTime_ms;
  }
}

//connect to the gyro sensor
//and then estimate the gyro offset value
void gyroSetup()
{
  // Initialize IMU.
  Wire.begin();
  if (!imu.init())
  {
    while(true)
    {
      Serial.println("Failed to detect and initialize IMU!");
      delay(200);
    }
  }
  imu.enableDefault();
  imu.writeReg(LSM6::CTRL2_G,  0b01011000); // 208 Hz, 1000 deg/s
  imu.writeReg(LSM6::CTRL1_XL, 0b01011000); // 208 Hz, +/- 4g

  // Wait for IMU readings to stabilize.
  delay(1000);

  // Calibrate the gyro (find the offset by taking an average of values)
  int32_t totalY = 0;
  for (int i = 0; i < CALIBRATION_ITERATIONS; i++)
  {
    imu.read();
    totalY += imu.g.y;
    delay(1);
  }

  gYZero = totalY / CALIBRATION_ITERATIONS;
  angleIntegrated = (int32_t) (atan2f((double) imu.a.z, (double) imu.a.x)*57295);

}

//imu.g.y is the rotation rate in mdps/LSB
//mdps: mili dgrees per second
//LSB: least significant bit of the data
//Depending on the setting, (imu.writeReg)
//the bit to mdps ratio will vary (see page 15 of LSM6DS33 manual)
//For settings chosen here, one increment corresponds to 35 mdps
//thus, to convert to degrees per second, multiply imu.g.y by 35/1000 ~ 1/29
void integrateGyro()
{
  // Convert from full-scale 1000 deg/s to deg/s.
  gY = (imu.g.y - gYZero) / 29;
  aX = imu.a.x;
  aZ = imu.a.z;
  // millidegrees
  //factor of 1000 is because dt is in ms
  angleIntegrated += gY * (curTime_ms - prevTime_ms); 
  angle_rad = ((float) angleIntegrated)/1000/180*3.14159;
  angleGravity = atan2f((double) aZ, (double) aX);

  angleGravityTemp = angleGravity;
  while(angleGravityTemp-angle_rad>=3.14159){
    angleGravityTemp-=6.28319;
  }

  while(angleGravityTemp-angle_rad<=-3.14159){
    angleGravityTemp+=6.28319;
  }

  angle_rad = (1-alpha)*angle_rad + alpha*angleGravityTemp;
}
