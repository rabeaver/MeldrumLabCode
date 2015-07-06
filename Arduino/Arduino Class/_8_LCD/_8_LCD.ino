#include <LiquidCrystal.h>
//ground 1,5,16
//power 5v 2,15
//pot 3
int i=1;
LiquidCrystal lcd(12, 11, 5, 4, 3, 2);
//                rs  en  d4 d5 d6 d7
void setup() {
  // set up the LCD's number of columns and rows: 
  lcd.begin(16, 2);
  // Print a message to the LCD.
  Serial.begin(9600);
}

void loop()
{
  lcd.setCursor(5,0);
  lcd.print("hello");
  lcd.setCursor(4,1);
  lcd.print(i);
  delay(500);
  i++;
}
