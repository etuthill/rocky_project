#include <Wire.h>
#include <LSM6.h>

LSM6 imu;

float angle_sum = 0;
int sample_count = 0;

void setup() {
  Serial.begin(9600);
  Wire.begin();

  if (!imu.init()) {
    Serial.println("IMU failed");
    while (1);
  }

  imu.enableDefault();

  Serial.println("Hold Rocky upright and still...");
  delay(3000);
}

void loop() {

  imu.read();

  // compute tilt using accelerometer
  float ax = imu.a.x;
  float az = imu.a.z;

  float angle = atan2(ax, az);  // radians

  Serial.print("angle sample: ");
  Serial.println(angle,6);

  angle_sum += angle;
  sample_count++;

  delay(10);

  if (sample_count >= 500) {

    float avg_angle = angle_sum / sample_count;

    Serial.println();
    Serial.println("===== RESULT =====");
    Serial.print("Average upright angle: ");
    Serial.println(avg_angle,6);

    Serial.println();
    Serial.print("#define FIXED_ANGLE_CORRECTION (");
    Serial.print(avg_angle,6);
    Serial.println(")");
    Serial.println("==================");

    while(1);
  }
}