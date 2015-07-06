//Brown = Ground
//Red = 5V power
//Yellow = Signal(Pin 9

#include <Servo.h>

Servo myservo;

int pos=0;

void setup() {
  // put your setup code here, to run once:
myservo.attach(9);

Serial.begin(9600);
}

void loop() {
 
  myservo.write(10);
  delay(500);
  
  myservo.write(10);
  delay(500);
  myservo.write(20);
  delay(500);
  myservo.write(30);
  delay(1000);
  myservo.write(40);
  delay(500);
  myservo.write(50);
  delay(500);
  myservo.write(60);
  delay(1000);
  myservo.write(70);
  delay(500);
  myservo.write(80);
  delay(500);

}

//Can you write your birthday in degrees? 
//Can you point to a semi-circle of 10 digits, and communicate your telephone number? 
