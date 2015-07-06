//Setup: Two pushbuttons, pull-HIGH resistors, 10k
//Connect high end of buttons to 7 + 8
//LED or Piezo(b) from pin 3

int toner = 500;


int Cscale[] = {262,294,330,349,392,440,494,523};
//               0   1   2   3   4   5   6   7
//      ~~       C   D   E   F   G   A   B   C
//call with Cscale[#] => freq
//not actually called in this script, but feel free to use!


void setup() {
  // put your setup code here, to run once:
pinMode(7, INPUT);
pinMode(8, INPUT);
pinMode(3, OUTPUT);
Serial.begin(9600);
}

void loop() {
  // put your main code here, to run repeatedly:
if(toner>2500){
  toner=300;}
if(toner<300){
  toner=2480;}

if(digitalRead(7) == LOW){
  toner +=20;
}
if(digitalRead(8) == LOW){
  toner -=20;
}
tone(3, toner);
//tone(3, Cscale[0=<i<8]);


Serial.println(toner);
delay(50);

}

//Can you play mary had a little lamb? Hint: A4 = 440hz
//Extra argument for tone(pin, freq, duration(millis))
//Come up with another input for this and map to some frequencies
