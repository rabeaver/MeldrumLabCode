//outputs go to the three shorter terminals of the RGB LED
//Blue is furthest from long stem, Green is second longest, Red is outside from long stem

void setup() {
pinMode(9,OUTPUT);  //blue
pinMode(10,OUTPUT); //green
pinMode(11,OUTPUT); //red
}



void loop() {

//we can light certain colors by calling something like
//analogWrite(9,255);
//analogWrite(10,0);
//analogWrite(11,255);
//it's purple! 
  
//but let's play a little bit:

for(int y=9; y<12;y++){

  for(int x=0; x<255; x++){
    analogWrite(y, x);
    delay(10); //fade from 0->255 for one color
  }
  analogWrite(y,0); //turn previous color off before moving on
}
}


//now make a LIGHTSHOW! 
