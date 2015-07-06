//Setup: Two pushbuttons, pull-HIGH resistors, 10k
//Connect high end of buttons to 7 + 8
//LED or Piezo(b) from pin 3

int brightness = 128;


void setup() {
  pinMode(7, INPUT);
  pinMode(8, INPUT);
  pinMode(3, OUTPUT);

  Serial.begin(9600);
}

void loop() {

if(digitalRead(7) == LOW){
  brightness +=2;
}
if(digitalRead(8) == LOW){
  brightness -=2;
}
constrain(brightness, 0, 255);

analogWrite(3, brightness);
Serial.println(brightness);
delay(5);

//try stopping it at one iteration per click
//hint: store the last state as an integer (0 or 1)
//and only take a change from low to high
}
