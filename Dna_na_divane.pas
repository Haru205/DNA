program divannoe_dna;

CONST AMINO_ACID = 'ARNDCQEGHILKMFPSTWYV';
//user crt;
type
    arr_t = array [1..3333] of char; //тип для последовательности аминокислот, в 
    //которую мы переводим заданный при вводе РНК
    arra_t = array [1..10000] of char; //тип для РНК
    dnk_t = array [1..500] of char; //тип для ДНК
    rubbish_t = array [1..10000] of integer;
    
    {Эта функция получает на вход кодон(последовательность из трех букв A,G,U,C,
    потом переводит эти три буквы в верхний регистр для удобства перевода, 
    каждому набору из трех букв ставит в соответствие некоторую другую букву, 
    обозначающую аминокислоту (по таблице)}
function Translate (codon:string) : char;
begin
    Translate := ' ';
    codon[1] := Upcase(codon[1]);
    codon[2] := Upcase(codon[2]);
    codon[3] := Upcase(codon[3]);
    if (codon = 'GCU') or (codon = 'GCC') or 
       (codon = 'GCA') or (codon = 'GCG') then
        Translate := 'A';
    if (codon = 'CGU') or (codon = 'CGC') or (codon = 'CGA') or
       (codon = 'CGG') or (codon = 'AGA') or (codon = 'AGG') then
        Translate := 'R';
    if (codon = 'AAU') or (codon = 'AAC') then
        Translate := 'N';
    if (codon = 'GAU') or (codon = 'GAC') then
        Translate := 'D';
    if (codon = 'UGU') or (codon = 'UGC') then
        Translate := 'C';
    if (codon = 'CAA') or (codon = 'CAG') then
        Translate := 'Q';
    if (codon = 'GAA') or (codon = 'GAG') then
        Translate := 'E';
    if (codon = 'GGU') or (codon = 'GGC') or 
       (codon = 'GGA') or (codon = 'GGG') then
        Translate := 'G';
    if (codon = 'CAU') or (codon = 'CAC') then
        Translate := 'H';
    if (codon = 'AUU') or (codon = 'AUC') or (codon = 'AUA') then
        Translate := 'I';
    if (codon = 'UUA') or (codon = 'UUG') or (codon = 'CUU') or
       (codon = 'CUC') or (codon = 'CUA') or (codon = 'CUG') then
        Translate := 'L';
    if (codon = 'AAA') or (codon = 'AAG') then
        Translate := 'K';
    if (codon = 'AUG') then
        Translate := 'M';
    if (codon = 'UUU') or (codon = 'UUC') then
        Translate := 'F';
    if (codon = 'CCU') or (codon = 'CCC') or 
       (codon = 'CCA') or (codon = 'CCG') then
        Translate := 'P';
    if (codon = 'UCU') or (codon = 'UCC') or (codon = 'UCA') or
       (codon = 'UCG') or (codon = 'AGU') or (codon = 'AGC') then
        Translate := 'S';
    if (codon = 'ACU') or (codon = 'ACC') or 
       (codon = 'ACA') or (codon = 'ACG') then
        Translate := 'T';
    if (codon = 'UGG') then
        Translate := 'W';
    if (codon = 'UAU') or (codon = 'UAC') then
        Translate := 'Y';
    if (codon = 'GUU') or (codon = 'GUC') or 
       (codon = 'GUA') or (codon = 'GUG') then
        Translate := 'V';
    {if (codon = 'UAA') or (codon = 'UGA') or
       (codon = 'UAG') then
         Translate := '0';}
end;

 {Эта функция принимает на вход последовательность РНК, а также длину РНК.
 После этого по три символа вводит в last_ch, передает ее в функцию Translate,
 чтобы сделать аминокислоту}
function Make_protein1 (rna:arra_t; index : integer) : arr_t;
var 
    index_trans : integer;
    last_ch : string;
    translation1 : arr_t;
begin
    index_trans:=1;
    while (index_trans <= index) do
    begin
		    last_ch := rna[index_trans];
		    if (index_trans + 1 <= index) then
		        last_ch := last_ch + rna[index_trans+1];
		    if (index_trans + 2 <= index) then
			      last_ch := last_ch + rna[index_trans+2];
        if (length(last_ch) = 3) then
			      translation1[(index_trans div 3)+1] := Translate(last_ch);
    index_trans := index_trans + 3;
	  end;
    Make_protein1 := translation1;   
end;

{Почему заданы три различные функции, которые делают почти одно и то же?
Это из-за того, что искомая аминокислота может начинаться на с первого символа
в РНК. Пример: пусть YYY-искомая последовательность. Может быть так: xxxYYY...,
xxYYY..., xYYY, то есть начало искомой последовательности дает разные остатки
по модулю 3}
function Make_protein2 (rna:arra_t; index : integer) : arr_t;
var 
    index_trans : integer;
    last_ch : string;
    translation2 : arr_t;
begin
    index_trans:=2;
    while (index_trans <= index) do
    begin
        last_ch := rna[index_trans];
        if (index_trans + 1 <= index) then
            last_ch := last_ch + rna[index_trans+1];
        if (index_trans + 2 <= index) then
            last_ch := last_ch + rna[index_trans+2];
        if (length(last_ch) = 3) then
            translation2[(index_trans div 3)+1] := Translate(last_ch);
        index_trans := index_trans + 3;
    end;
    Make_protein2 := translation2;   
end;    

function Make_protein3 (rna:arra_t; index : integer) : arr_t;
var 
    index_trans : integer;
    last_ch : string;
    translation3 : arr_t;
begin
    index_trans:=3;
    while (index_trans <= index) do
    begin
        last_ch := rna[index_trans];
        if (index_trans + 1 <= index) then
            last_ch := last_ch + rna[index_trans+1];
        if (index_trans + 2 <= index) then
            last_ch := last_ch + rna[index_trans+2];
        if (length(last_ch) = 3) then
            translation3[index_trans div 3] := Translate(last_ch);
        index_trans := index_trans + 3;
    end;
    Make_protein3 := translation3;   
end;

{Небольшая функция, в которую подается три параметра: индекс первого вхождения
заданного в условии белка в последовательности РНК из условия со сдвигом 
(см.предыдущие функции)}
function Sortirovka(znach1,znach2,znach3:integer) : integer;
var 
    temp1:integer;
begin
    if(znach1>znach2) then
    begin
        temp1 := znach1;
        znach1 := znach2;
        znach2 := temp1;
    end;
    if(znach2>znach3) then
    begin
        temp1 := znach2;
        znach2 := znach3;
        znach3 := temp1;
    end;
    if(znach1>znach2) then
    begin
        temp1 := znach1;
        znach1 := znach2;
        znach2 := temp1;
    end;
    if(znach1<>-1) then
        Sortirovka := znach1
    else 
        if (znach2<>-1) then
            Sortirovka := znach2
        else
            if(znach3<>-1) then
                Sortirovka := znach3
            else 
                Sortirovka := -1;
end;

{Функция возвращает индекс первого вхождения белка в последовательность 
аминокислот, которую мы получили из РНК, иначе -1}
function Number(trans:arr_t; dna: dnk_t;index,index_dna: integer) : integer;
var 
    num, ind,temp : integer;
begin
    num := -1; temp:=1;
    while (num=(-1)) and (temp<=index-index_dna+1) do
    begin
        num:=temp;
        for ind:=1 to index_dna do
            if(trans[temp+ind-1] <> dna [ind]) then
                num := -1;
        if(num<>-1) then
            break;
            temp := temp+1;
    end;
    Number := num;
end;

{Главная функция, в ней вводятся последовательность РНК и белок, обрабаотываются,
с помощью дополнительных функций ищет индекс первого вхождения заданного белка
в последовательность РНК}
procedure Main ();
var
    translation1,translation2,translation3 : arr_t; {Три различных последовательности
    Аминокислот из РНК}
    new_character : char; {Очередной вводимый символ}
    znach1, znach2, znach3, otvet : integer; {Три индекса первого вхождения
    в рахных последовательностях и ответ}
    index1,index2,index3 : integer;
    rna : arra_t; // РНК
    rubbish_count : rubbish_t; {Массив, в котором считается количество
    пробелов и тире между символами}
    dna : dnk_t; // Белок
    index, index_dna : integer; {Длина РНК и длина белка}
    error : string; {Ловец ошибок}
    nachalo, konec, count : integer; {Вывод начала и конца в ответ, счетчик}
begin 
    writeln('Vvedit RNA and Belok cherez tochky');
    read(new_character);
    count := 0;
    index := 0;
    znach1 := -1;znach2 := -1;znach3 := -1;
    error := '';
    index_dna := 0;
    {Ввод последовательности РНК}
	  while (new_character <> '.') and (index <= 500) and (error <> 'error') do
	      begin
	      if (index<=10000) then
	      begin
			      if(new_character = 'A') or (new_character = 'a') or
	            (new_character = 'U') or (new_character = 'u') or
	            (new_character = 'G') or (new_character = 'g') or
	            (new_character = 'C') or (new_character = 'c') then
	          begin
	              index := index + 1;
	              rubbish_count[index] := count;
	              count := 0;
	              rna[index] := new_character;
	          end
		        else
	          begin
		            if (new_character = '-') or (new_character = ' ') or
		               (new_character = '.') then
                begin
                    if (new_character = '-') or (new_character = ' ') then
			                    count := count + 1
                end
                else
			            error := 'error';  
		        end;
		    end
		    else
		        error := 'error';
	      read(new_character);
	      end;
	      
	      {Ввод белка}
	  if(error <> 'error')  then
    begin    
        read(new_character);
	      while (new_character <> '.')  and (error <> 'error')and (new_character<>'') do
	      begin    
	          if (index_dna<=500) then
	          begin
			          if(Pos(Upcase(new_character),amino_acid) <> 0) then
	              begin
	                  index_dna := index_dna + 1;
	                  dna[index_dna] := Upcase(new_character);
	              end
		            else
	              begin
		                if (new_character <> '-') and (new_character <> ' ') then
			                  error := 'error';  
		            end;
		        end
		        else
		            error := 'error';
        read(new_character);
        end;
    end;
    {Обработка РНК, перевод в последовательности аминокислот}
    if (new_character = ' ') and (error<>'error')or (index_dna = 0) then
      writeln('Incorrect input')
    else
    begin
        translation1 := Make_protein1(rna,index);
	      translation2 := Make_protein2(rna,index);
	      translation3 := Make_protein3(rna,index);
	      if ((index mod 3) = 0) then
        begin
            index1 := index div 3;
            index2 := index div 3 - 1;
            index3 := index div 3 - 1;
        end;
        if ((index mod 3) = 1) then
        begin
            index1 := index div 3;
            index2 := index div 3;
            index3 := index div 3 - 1;
        end;
        if ((index mod 3) = 2) then
        begin
            index1 := index div 3;
            index2 := index div 3;
            index3 := index div 3;
        end;
        if (index_dna>index div 3) then
            writeln('Nichego ne naideno')
        else
          {Поиск ответа, затем его вывод}
        begin
            znach1 := Number(translation1,dna,index1,index_dna);
            if (znach1<>-1) then
                znach1 := 3*(znach1-1)+1;
            znach2 := Number(translation2,dna,index2,index_dna);
            if (znach2<>-1) then
                znach2 := (znach2-1)*3+2;
            znach3 := Number(translation3,dna,index3,index_dna);
            if (znach3<>-1) then
                znach3 := znach3*3;
            otvet := Sortirovka(znach1,znach2,znach3);
            if (otvet = -1) then
                writeln('Nichego ne naideno');
            if (otvet <> -1) then
            begin
                nachalo := 0;
                konec := 0;
                for index1:= 1 to otvet do
                    nachalo := nachalo + rubbish_count[index1];
                for index1:= 1 to otvet+3*index_dna-1 do
                    konec := konec + rubbish_count[index1];
                writeln('V posledovatelnosti so znakami plus i tire');
                writeln('Start: ',nachalo+otvet,'   End: ', konec+otvet+3*index_dna-1);
                writeln('In posledovatelnostt:');
                writeln('Start: ',otvet,'   End: ',otvet+3*index_dna-1);
                writeln('Iskomaya posledovatelnostt:');
                if (znach1 = otvet) then
                    for znach1:= otvet to (otvet+3*index_dna-1) do
                        write(rna[znach1]);
                if (znach2 = otvet) then
                    for znach2:= otvet to (otvet+3*index_dna-1) do
                        write(rna[znach2]);
                if (znach3 = otvet) then
                    for znach3:= otvet to (otvet+3*index_dna-1) do
                        write(rna[znach3]);
                writeln();
            end;  
        end;
    end;
end;

var
    iteration : char;
    stroka, oshibka : string;
    schetcik, condition : integer;
begin
    condition := 1;
    {Повторять ввода РНК и белка много раз, пока последний ввод не будет равен
    end}
    while (condition = 1) do
    begin
        Main();
        writeln('Vvedite "end."(with tochka), esli want zakonchit');
        writeln('Inache vvedite tochky');
        read(iteration);
        stroka := '';
        while (iteration<>'.') do
        begin
            stroka := stroka+iteration;
            read(iteration);
            end;
        if(Pos('end',stroka)>0) then
            condition := 0;
    end;
end.