/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...

#include "Riostream.h"

#include "fitF.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"

ClassImp(fitF)

 fitF::fitF(const char *name, const char *title,
                        RooAbsReal& _x,
                        RooAbsCategory& _y,
                        RooAbsReal& _neueff,
                        RooRealVar *_pp[]
            ) :
   RooAbsPdf(name,title),
   x("x","x",this,_x),
   y("y","y",this,_y),
   neueff("neueff","neueff",this,_neueff),
   p0("p0","p0",this,*_pp[0]),//l0
   p1("p1","p1",this,*_pp[1]),
   p2("p2","p2",this,*_pp[2]),
   p3("p3","p3",this,*_pp[3]),
   p4("p4","p4",this,*_pp[4]),
   p5("p5","p5",this,*_pp[5]),
   p6("p6","p6",this,*_pp[6]),
   p7("p7","p7",this,*_pp[7]),
   p8("p8","p8",this,*_pp[8]),
   p9("p9","p9",this,*_pp[9]),//p1n
   p10("p10","p10",this,*_pp[10]),
   p11("p11","p11",this,*_pp[11]),
   p12("p12","p12",this,*_pp[12]),
   p13("p13","p13",this,*_pp[13]),
   p14("p14","p14",this,*_pp[14]),
   p15("p15","p15",this,*_pp[15]),
   p16("p16","p16",this,*_pp[16]),
   p17("p17","p17",this,*_pp[17]),
   p18("p18","p18",this,*_pp[18]),//p2n
   p19("p19","p19",this,*_pp[19]),
   p20("p20","p20",this,*_pp[20]),
   p21("p21","p21",this,*_pp[21]),
   p22("p22","p22",this,*_pp[22]),
   p23("p23","p23",this,*_pp[23]),
   p24("p24","p24",this,*_pp[24]),
   p25("p25","p25",this,*_pp[25]),
   p26("p26","p26",this,*_pp[26]),
   p27("p27","p27",this,*_pp[27]),
   p28("p28","p28",this,*_pp[28]),
   p29("p29","p29",this,*_pp[29]),
   p30("p30","p30",this,*_pp[30]),
   p31("p31","p31",this,*_pp[31]),
   p32("p32","p32",this,*_pp[32]),
   p33("p33","p33",this,*_pp[33]),
   p34("p34","p34",this,*_pp[34]),
   p35("p35","p35",this,*_pp[35]),
   p36("p36","p36",this,*_pp[36]),
   p37("p37","p37",this,*_pp[37]),
   p38("p38","p38",this,*_pp[38]),
   p39("p39","p39",this,*_pp[39]),
   p40("p40","p40",this,*_pp[40]),
   p41("p41","p41",this,*_pp[41]),
   p42("p42","p42",this,*_pp[42]),
   p43("p43","p43",this,*_pp[43]),
   p44("p44","p44",this,*_pp[44]),
   p45("p45","p45",this,*_pp[45]),
   p46("p46","p46",this,*_pp[46]),
   p47("p47","p47",this,*_pp[47]),
   p48("p48","p48",this,*_pp[48]),
   p49("p49","p49",this,*_pp[49]),
   p50("p50","p50",this,*_pp[50]),
   p51("p51","p51",this,*_pp[51]),
   p52("p52","p52",this,*_pp[52]),
   p53("p53","p53",this,*_pp[53]),
   p54("p54","p54",this,*_pp[54]),
   p55("p55","p55",this,*_pp[55]),
   p56("p56","p56",this,*_pp[56]),
   p57("p57","p57",this,*_pp[57]),
   p58("p58","p58",this,*_pp[58]),
   p59("p59","p59",this,*_pp[59]),
   p60("p60","p60",this,*_pp[60]),
   p61("p61","p61",this,*_pp[61]),
   p62("p62","p62",this,*_pp[62]),
   p63("p63","p63",this,*_pp[63]),
   p64("p64","p64",this,*_pp[64]),
   p65("p65","p65",this,*_pp[65]),
   p66("p66","p66",this,*_pp[66]),
   p67("p67","p67",this,*_pp[67]),
   p68("p68","p68",this,*_pp[68]),
   p69("p69","p69",this,*_pp[69]),
   p70("p70","p70",this,*_pp[70]),
   p71("p71","p71",this,*_pp[71]),
   p72("p72","p72",this,*_pp[72]),
   p73("p73","p73",this,*_pp[73]),
   p74("p74","p74",this,*_pp[74]),
   p75("p75","p75",this,*_pp[75]),
   p76("p76","p76",this,*_pp[76]),
   p77("p77","p77",this,*_pp[77]),
   p78("p78","p78",this,*_pp[78]),
   p79("p79","p79",this,*_pp[79]),
   p80("p80","p80",this,*_pp[80]),
   p81("p81","p81",this,*_pp[81]),
   p82("p82","p82",this,*_pp[82]),
   p83("p83","p83",this,*_pp[83]),
   p84("p84","p84",this,*_pp[84]),
   p85("p85","p85",this,*_pp[85]),
   p86("p86","p86",this,*_pp[86]),
   p87("p87","p87",this,*_pp[87]),
   p88("p88","p88",this,*_pp[88]),
   p89("p89","p89",this,*_pp[89]),
   p90("p90","p90",this,*_pp[90]),
   p91("p91","p91",this,*_pp[91]),
   p92("p92","p92",this,*_pp[92]),
   p93("p93","p93",this,*_pp[93]),
   p94("p94","p94",this,*_pp[94]),
   p95("p95","p95",this,*_pp[95]),
   p96("p96","p96",this,*_pp[96]),
   p97("p97","p97",this,*_pp[97]),
   p98("p98","p98",this,*_pp[98]),
   p99("p99","p99",this,*_pp[99]),
   p100("p100","p100",this,*_pp[100]),
   p101("p101","p101",this,*_pp[101]),
   p102("p102","p102",this,*_pp[102]),
   p103("p103","p103",this,*_pp[103]),
   p104("p104","p104",this,*_pp[104]),
   p105("p105","p105",this,*_pp[105]),

   p106("p106","p106",this,*_pp[106]),
   p107("p107","p107",this,*_pp[107]),
   p108("p108","p108",this,*_pp[108]),
   p109("p109","p109",this,*_pp[109]),
   p110("p110","p110",this,*_pp[110]),
   p111("p111","p111",this,*_pp[111]),
   p112("p112","p112",this,*_pp[112]),
   p113("p113","p113",this,*_pp[113]),
   p114("p114","p114",this,*_pp[114]),
   p115("p115","p115",this,*_pp[115]),
   p116("p116","p116",this,*_pp[116]),
   p117("p117","p117",this,*_pp[117]),
   p118("p118","p118",this,*_pp[118]),
   p119("p119","p119",this,*_pp[119]),
   p120("p120","p120",this,*_pp[120]),
   p121("p121","p121",this,*_pp[121]),
   p122("p122","p122",this,*_pp[122]),
   p123("p123","p123",this,*_pp[123]),
   p124("p124","p124",this,*_pp[124]),
   p125("p125","p125",this,*_pp[125]),
   p126("p126","p126",this,*_pp[126]),
   p127("p127","p127",this,*_pp[127]),
   p128("p128","p128",this,*_pp[128]),
   p129("p129","p129",this,*_pp[129]),
   p130("p130","p130",this,*_pp[130]),
   p131("p131","p131",this,*_pp[131]),
   p132("p132","p132",this,*_pp[132]),
   p133("p133","p133",this,*_pp[133]),
   p134("p134","p134",this,*_pp[134]),
   p135("p135","p135",this,*_pp[135]),
   p136("p136","p136",this,*_pp[136]),
   p137("p137","p137",this,*_pp[137]),
   p138("p138","p138",this,*_pp[138]),
   p139("p139","p139",this,*_pp[139]),
   p140("p140","p140",this,*_pp[140]),
   p141("p141","p141",this,*_pp[141]),
   p142("p142","p142",this,*_pp[142]),
   p143("p143","p143",this,*_pp[143]),
   p144("p144","p144",this,*_pp[144]),
   p145("p145","p145",this,*_pp[145]),
   p146("p146","p146",this,*_pp[146]),
   p147("p147","p147",this,*_pp[147]),
   p148("p148","p148",this,*_pp[148]),
   p149("p149","p149",this,*_pp[149]),
   p150("p150","p150",this,*_pp[150])
 {


 }


 fitF::fitF(const fitF& other, const char* name) :
   RooAbsPdf(other,name),
   x("x",this,other.x),
   y("y",this,other.y),
   neueff("neueff",this,other.neueff),
   p0("p0",this,other.p0),
   p1("p1",this,other.p1),
   p2("p2",this,other.p2),
   p3("p3",this,other.p3),
   p4("p4",this,other.p4),
   p5("p5",this,other.p5),
   p6("p6",this,other.p6),
   p7("p7",this,other.p7),
   p8("p8",this,other.p8),
   p9("p9",this,other.p9),
   p10("p10",this,other.p10),
   p11("p11",this,other.p11),
   p12("p12",this,other.p12),
   p13("p13",this,other.p13),
   p14("p14",this,other.p14),
   p15("p15",this,other.p15),
   p16("p16",this,other.p16),
   p17("p17",this,other.p17),
   p18("p18",this,other.p18),
   p19("p19",this,other.p19),
   p20("p20",this,other.p20),
   p21("p21",this,other.p21),
   p22("p22",this,other.p22),
   p23("p23",this,other.p23),
   p24("p24",this,other.p24),
   p25("p25",this,other.p25),
   p26("p26",this,other.p26),
   p27("p27",this,other.p27),
   p28("p28",this,other.p28),
   p29("p29",this,other.p29),
   p30("p30",this,other.p30),
   p31("p31",this,other.p31),
   p32("p32",this,other.p32),
   p33("p33",this,other.p33),
   p34("p34",this,other.p34),
   p35("p35",this,other.p35),
   p36("p36",this,other.p36),
   p37("p37",this,other.p37),
   p38("p38",this,other.p38),
   p39("p39",this,other.p39),
   p40("p40",this,other.p40),
   p41("p41",this,other.p41),
   p42("p42",this,other.p42),
   p43("p43",this,other.p43),
   p44("p44",this,other.p44),
   p45("p45",this,other.p45),
   p46("p46",this,other.p46),
   p47("p47",this,other.p47),
   p48("p48",this,other.p48),
   p49("p49",this,other.p49),
   p50("p50",this,other.p50),
   p51("p51",this,other.p51),
   p52("p52",this,other.p52),
   p53("p53",this,other.p53),
   p54("p54",this,other.p54),
   p55("p55",this,other.p55),
   p56("p56",this,other.p56),
   p57("p57",this,other.p57),
   p58("p58",this,other.p58),
   p59("p59",this,other.p59),
   p60("p60",this,other.p60),
   p61("p61",this,other.p61),
   p62("p62",this,other.p62),
   p63("p63",this,other.p63),
   p64("p64",this,other.p64),
   p65("p65",this,other.p65),
   p66("p66",this,other.p66),
   p67("p67",this,other.p67),
   p68("p68",this,other.p68),
   p69("p69",this,other.p69),
   p70("p70",this,other.p70),
   p71("p71",this,other.p71),
   p72("p72",this,other.p72),
   p73("p73",this,other.p73),
   p74("p74",this,other.p74),
   p75("p75",this,other.p75),
   p76("p76",this,other.p76),
   p77("p77",this,other.p77),
   p78("p78",this,other.p78),
   p79("p79",this,other.p79),
   p80("p80",this,other.p80),
   p81("p81",this,other.p81),
   p82("p82",this,other.p82),
   p83("p83",this,other.p83),
   p84("p84",this,other.p84),
   p85("p85",this,other.p85),
   p86("p86",this,other.p86),
   p87("p87",this,other.p87),
   p88("p88",this,other.p88),
   p89("p89",this,other.p89),
   p90("p90",this,other.p90),
   p91("p91",this,other.p91),
   p92("p92",this,other.p92),
   p93("p93",this,other.p93),
   p94("p94",this,other.p94),
   p95("p95",this,other.p95),
   p96("p96",this,other.p96),
   p97("p97",this,other.p97),
   p98("p98",this,other.p98),
   p99("p99",this,other.p99),
   p100("p100",this,other.p100),
   p101("p101",this,other.p101),
   p102("p102",this,other.p102),
   p103("p103",this,other.p103),
   p104("p104",this,other.p104),
   p105("p105",this,other.p105),

   p106("p106",this,other.p106),
   p107("p107",this,other.p107),
   p108("p108",this,other.p108),
   p109("p109",this,other.p109),
   p110("p110",this,other.p110),
   p111("p111",this,other.p111),
   p112("p112",this,other.p112),
   p113("p113",this,other.p113),
   p114("p114",this,other.p114),
   p115("p115",this,other.p115),
   p116("p116",this,other.p116),
   p117("p117",this,other.p117),
   p118("p118",this,other.p118),
   p119("p119",this,other.p119),
   p120("p120",this,other.p120),
   p121("p121",this,other.p121),
   p122("p122",this,other.p122),
   p123("p123",this,other.p123),
   p124("p124",this,other.p124),
   p125("p125",this,other.p125),
   p126("p126",this,other.p126),
   p127("p127",this,other.p127),
   p128("p128",this,other.p128),
   p129("p129",this,other.p129),
   p130("p130",this,other.p130),
   p131("p131",this,other.p131),
   p132("p132",this,other.p132),
   p133("p133",this,other.p133),
   p134("p134",this,other.p134),
   p135("p135",this,other.p135),
   p136("p136",this,other.p136),
   p137("p137",this,other.p137),
   p138("p138",this,other.p138),
   p139("p139",this,other.p139),
   p140("p140",this,other.p140),
   p141("p141",this,other.p141),
   p142("p142",this,other.p142),
   p143("p143",this,other.p143),
   p144("p144",this,other.p144),
   p145("p145",this,other.p145),
   p146("p146",this,other.p146),
   p147("p147",this,other.p147),
   p148("p148",this,other.p148),
   p149("p149",this,other.p149),
   p150("p150",this,other.p150)
 {

 }



 Double_t fitF::evaluate() const
 {

     Double_t par0[151];
     par0[0]=p0;
     par0[1]=p1;
     par0[2]=p2;
     par0[3]=p3;
     par0[4]=p4;
     par0[5]=p5;
     par0[6]=p6;
     par0[7]=p7;
     par0[8]=p8;
     par0[9]=p9;
     par0[10]=p10;
     par0[11]=p11;
     par0[12]=p12;
     par0[13]=p13;
     par0[14]=p14;
     par0[15]=p15;
     par0[16]=p16;
     par0[17]=p17;
     par0[18]=p18;

     par0[19]=p19;
     par0[20]=p20;
     par0[21]=p21;
     par0[22]=p22;
     par0[23]=p23;
     par0[24]=p24;
     par0[25]=p25;
     par0[26]=p26;
     par0[27]=p27;
     par0[28]=p28;
     par0[29]=p29;
     par0[30]=p30;
     par0[31]=p31;
     par0[32]=p32;
     par0[33]=p33;
     par0[34]=p34;
     par0[35]=p35;
     par0[36]=p36;
     par0[37]=p37;
     par0[38]=p38;
     par0[39]=p39;
     par0[40]=p40;
     par0[41]=p41;
     par0[42]=p42;
     par0[43]=p43;
     par0[44]=p44;
     par0[45]=p45;
     par0[46]=p46;
     par0[47]=p47;
     par0[48]=p48;
     par0[49]=p49;
     par0[50]=p50;
     par0[51]=p51;
     par0[52]=p52;
     par0[53]=p53;
     par0[54]=p54;
     par0[55]=p55;
     par0[56]=p56;
     par0[57]=p57;
     par0[58]=p58;
     par0[59]=p59;
     par0[60]=p60;
     par0[61]=p61;
     par0[62]=p62;
     par0[63]=p63;
     par0[64]=p64;
     par0[65]=p65;
     par0[66]=p66;
     par0[67]=p67;
     par0[68]=p68;
     par0[69]=p69;
     par0[70]=p70;
     par0[71]=p71;
     par0[72]=p72;
     par0[73]=p73;
     par0[74]=p74;
     par0[75]=p75;
     par0[76]=p76;
     par0[77]=p77;
     par0[78]=p78;
     par0[79]=p79;
     par0[80]=p80;
     par0[81]=p81;
     par0[82]=p82;
     par0[83]=p83;
     par0[84]=p84;
     par0[85]=p85;
     par0[86]=p86;
     par0[87]=p87;
     par0[88]=p88;
     par0[89]=p89;
     par0[90]=p90;
     par0[91]=p91;
     par0[92]=p92;
     par0[93]=p93;
     par0[94]=p94;
     par0[95]=p95;
     par0[96]=p96;
     par0[97]=p97;
     par0[98]=p98;
     par0[99]=p99;
     par0[100]=p100;
     par0[101]=p101;
     par0[102]=p102;
     par0[103]=p103;
     par0[104]=p104;
     par0[105]=p105;

     par0[106]=p106;
     par0[107]=p107;
     par0[108]=p108;
     par0[109]=p109;
     par0[110]=p110;
     par0[111]=p111;
     par0[112]=p112;
     par0[113]=p113;
     par0[114]=p114;
     par0[115]=p115;
     par0[116]=p116;
     par0[117]=p117;
     par0[118]=p118;
     par0[119]=p119;
     par0[120]=p120;
     par0[121]=p121;
     par0[122]=p122;
     par0[123]=p123;
     par0[124]=p124;
     par0[125]=p125;
     par0[126]=p126;
     par0[127]=p127;
     par0[128]=p128;
     par0[129]=p129;
     par0[130]=p130;
     par0[131]=p131;
     par0[132]=p132;
     par0[133]=p133;
     par0[134]=p134;
     par0[135]=p135;
     par0[136]=p136;
     par0[137]=p137;
     par0[138]=p138;
     par0[139]=p139;
     par0[140]=p140;
     par0[141]=p141;
     par0[142]=p142;
     par0[143]=p143;
     par0[144]=p144;
     par0[145]=p145;
     par0[146]=p146;
     par0[147]=p147;
     par0[148]=p148;
     par0[149]=p149;
     par0[150]=p150;


     Double_t par1[151];
     par1[0]=p0;
     par1[1]=p1;
     par1[2]=p2;
     par1[3]=p3;
     par1[4]=p4;
     par1[5]=p5;
     par1[6]=p6;
     par1[7]=p7;
     par1[8]=p8;
     par1[9]=p9;
     par1[10]=p10;
     par1[11]=p11;
     par1[12]=p12;
     par1[13]=p13;
     par1[14]=p14;
     par1[15]=p15;
     par1[16]=p16;
     par1[17]=p17;
     par1[18]=p18;
     par1[19]=p19;
     par1[20]=p20;

     par1[21]=p21;
     par1[22]=p22;
     par1[23]=p23;
     par1[24]=p24;
     par1[25]=p25;
     par1[26]=p26;
     par1[27]=p27;
     par1[28]=p28;
     par1[29]=p29;
     par1[30]=p30;
     par1[31]=p31;
     par1[32]=p32;
     par1[33]=p33;
     par1[34]=p34;
     par1[35]=p35;
     par1[36]=p36;
     par1[37]=p37;
     par1[38]=p38;
     par1[39]=p39;
     par1[40]=p40;
     par1[41]=p41;
     par1[42]=p42;
     par1[43]=p43;
     par1[44]=p44;
     par1[45]=p45;
     par1[46]=p46;
     par1[47]=p47;
     par1[48]=p48;
     par1[49]=p49;
     par1[50]=p50;
     par1[51]=p51;
     par1[52]=p52;
     par1[53]=p53;
     par1[54]=p54;
     par1[55]=p55;
     par1[56]=p56;
     par1[57]=p57;
     par1[58]=p58;
     par1[59]=p59;
     par1[60]=p60;
     par1[61]=p61;
     par1[62]=p62;
     par1[63]=p63;
     par1[64]=p64;
     par1[65]=p65;
     par1[66]=p66;
     par1[67]=p67;
     par1[68]=p68;
     par1[69]=p69;
     par1[70]=p70;
     par1[71]=p71;
     par1[72]=p72;
     par1[73]=p73;
     par1[74]=p74;
     par1[75]=p75;
     par1[76]=p76;
     par1[77]=p77;
     par1[78]=p78;
     par1[79]=p79;
     par1[80]=p80;
     par1[81]=p81;
     par1[82]=p82;
     par1[83]=p83;
     par1[84]=p84;
     par1[85]=p85;
     par1[86]=p86;
     par1[87]=p87;
     par1[88]=p88;
     par1[89]=p89;
     par1[90]=p90;
     par1[91]=p91;
     par1[92]=p92;
     par1[93]=p93;
     par1[94]=p94;
     par1[95]=p95;
     par1[96]=p96;
     par1[97]=p97;
     par1[98]=p98;
     par1[99]=p99;
     par1[100]=p100;
     par1[101]=p101;
     par1[102]=p102;
     par1[103]=p103;
     par1[104]=p104;
     par1[105]=p105;

     par1[106]=p106;
     par1[107]=p107;
     par1[108]=p108;
     par1[109]=p109;
     par1[110]=p110;
     par1[111]=p111;
     par1[112]=p112;
     par1[113]=p113;
     par1[114]=p114;
     par1[115]=p115;
     par1[116]=p116;
     par1[117]=p117;
     par1[118]=p118;
     par1[119]=p119;
     par1[120]=p120;
     par1[121]=p121;
     par1[122]=p122;
     par1[123]=p123;
     par1[124]=p124;
     par1[125]=p125;
     par1[126]=p126;
     par1[127]=p127;
     par1[128]=p128;
     par1[129]=p129;
     par1[130]=p130;
     par1[131]=p131;
     par1[132]=p132;
     par1[133]=p133;
     par1[134]=p134;
     par1[135]=p135;
     par1[136]=p136;
     par1[137]=p137;
     par1[138]=p138;
     par1[139]=p139;
     par1[140]=p140;
     par1[141]=p141;
     par1[142]=p142;
     par1[143]=p143;
     par1[144]=p144;
     par1[145]=p145;
     par1[146]=p146;
     par1[147]=p147;
     par1[148]=p148;
     par1[149]=p149;
     par1[150]=p150;

     Double_t par2[151];
     par2[0]=p0;
     par2[1]=p1;
     par2[2]=p2;
     par2[3]=p3;
     par2[4]=p4;
     par2[5]=p5;
     par2[6]=p6;
     par2[7]=p7;
     par2[8]=p8;
     par2[9]=p9;
     par2[10]=p10;
     par2[11]=p11;
     par2[12]=p12;
     par2[13]=p13;
     par2[14]=p14;
     par2[15]=p15;
     par2[16]=p16;
     par2[17]=p17;
     par2[18]=p18;
     par2[19]=p19;
     par2[20]=p20;
     par2[21]=p21;

     par2[22]=p22;
     par2[23]=p23;
     par2[24]=p24;
     par2[25]=p25;
     par2[26]=p26;
     par2[27]=p27;
     par2[28]=p28;
     par2[29]=p29;
     par2[30]=p30;
     par2[31]=p31;
     par2[32]=p32;
     par2[33]=p33;
     par2[34]=p34;
     par2[35]=p35;
     par2[36]=p36;
     par2[37]=p37;
     par2[38]=p38;
     par2[39]=p39;
     par2[40]=p40;
     par2[41]=p41;
     par2[42]=p42;
     par2[43]=p43;
     par2[44]=p44;
     par2[45]=p45;
     par2[46]=p46;
     par2[47]=p47;
     par2[48]=p48;
     par2[49]=p49;
     par2[50]=p50;
     par2[51]=p51;
     par2[52]=p52;
     par2[53]=p53;
     par2[54]=p54;
     par2[55]=p55;
     par2[56]=p56;
     par2[57]=p57;
     par2[58]=p58;
     par2[59]=p59;
     par2[60]=p60;
     par2[61]=p61;
     par2[62]=p62;
     par2[63]=p63;
     par2[64]=p64;
     par2[65]=p65;
     par2[66]=p66;
     par2[67]=p67;
     par2[68]=p68;
     par2[69]=p69;
     par2[70]=p70;
     par2[71]=p71;
     par2[72]=p72;
     par2[73]=p73;
     par2[74]=p74;
     par2[75]=p75;
     par2[76]=p76;
     par2[77]=p77;
     par2[78]=p78;
     par2[79]=p79;
     par2[80]=p80;
     par2[81]=p81;
     par2[82]=p82;
     par2[83]=p83;
     par2[84]=p84;
     par2[85]=p85;
     par2[86]=p86;
     par2[87]=p87;
     par2[88]=p88;
     par2[89]=p89;
     par2[90]=p90;
     par2[91]=p91;
     par2[92]=p92;
     par2[93]=p93;
     par2[94]=p94;
     par2[95]=p95;
     par2[96]=p96;
     par2[97]=p97;
     par2[98]=p98;
     par2[99]=p99;
     par2[100]=p100;
     par2[101]=p101;
     par2[102]=p102;
     par2[103]=p103;
     par2[104]=p104;
     par2[105]=p105;

     par2[106]=p106;
     par2[107]=p107;
     par2[108]=p108;
     par2[109]=p109;
     par2[110]=p110;
     par2[111]=p111;
     par2[112]=p112;
     par2[113]=p113;
     par2[114]=p114;
     par2[115]=p115;
     par2[116]=p116;
     par2[117]=p117;
     par2[118]=p118;
     par2[119]=p119;
     par2[120]=p120;
     par2[121]=p121;
     par2[122]=p122;
     par2[123]=p123;
     par2[124]=p124;
     par2[125]=p125;
     par2[126]=p126;
     par2[127]=p127;
     par2[128]=p128;
     par2[129]=p129;
     par2[130]=p130;
     par2[131]=p131;
     par2[132]=p132;
     par2[133]=p133;
     par2[134]=p134;
     par2[135]=p135;
     par2[136]=p136;
     par2[137]=p137;
     par2[138]=p138;
     par2[139]=p139;
     par2[140]=p140;
     par2[141]=p141;
     par2[142]=p142;
     par2[143]=p143;
     par2[144]=p144;
     par2[145]=p145;
     par2[146]=p146;
     par2[147]=p147;
     par2[148]=p148;
     par2[149]=p149;
     par2[150]=p150;

     Double_t t[1];
     t[0]=x;

     Double_t ret=0;
     if (y==0){
         ret = fcn_decay(t,par0)-fcn_1ndecay(t,par1)-fcn_2ndecay(t,par2);
     }else if (y==1){
         ret = fcn_1ndecay(t,par1);
     }else{
         ret = fcn_2ndecay(t,par2);
     }

     return ret;
 }

 //! Global function
 Double_t fitF::fcn_decay(Double_t *x, Double_t *par) const{
     Double_t* pn=&par[knri];
     Double_t* lamda=par;
     Double_t* p2n=&par[knri*2];
     Double_t* production_yield=&par[knri*3];

     //init
     Double_t N0=par[knri*4]/par[0];

     //! Parent nuclei
     Double_t returnval=lamda[0]*N0*exp(-lamda[0]*x[0]);

     for (Int_t i=0;i<npaths;i++){
         returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0]);
     }
     return returnval;
 }

 //! Function for decay with 1 neutron emission
 Double_t fitF::fcn_1ndecay(Double_t *x, Double_t *par) const{
     Double_t* pn=&par[knri];
     Double_t* lamda=par;
     Double_t* p2n=&par[knri*2];
     Double_t* production_yield=&par[knri*3];
     Double_t randcoinfgt0n=par[knri*4+2];
     Double_t randcoinf1n=par[knri*4+1];
     //init
     Double_t N0=par[knri*4]/par[0];

     //! Parent nuclei

     //! random coinc of beta decay of parent
     Double_t returnval=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;
     //! decay with 1 neutron of parent
     returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf1n-randcoinfgt0n);

     //! decay with 1 neutron of parent from p2n
     returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf1n-randcoinfgt0n);
     //! decay with 2 neutron of parent (not random 1 neutron)
     returnval-=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;

     for (Int_t i=0;i<npaths;i++){
         //! random coinc of beta decay of daugter
         returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*randcoinf1n;
         //! decay with 1 neutron of daugter nuclei
         returnval+=neueff*pn[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);

         //! decay with 1 neutron of daugter from p2n
         returnval+=2*(neueff*(1-neueff))*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);
         //! decay with 2 neutron of parent (not random 1 neutron)
         returnval-=neueff*neueff*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*randcoinf1n;
     }

     return returnval;

 }


 //! Function for decay with 2 neutron emission
 Double_t fitF::fcn_2ndecay(Double_t *x, Double_t *par) const{
     Double_t* pn=&par[knri];
     Double_t* lamda=par;
     Double_t* p2n=&par[knri*2];
     Double_t* production_yield=&par[knri*3];

     Double_t randcoinf2n=par[knri*4+3];
     Double_t randcoinfgt0n=par[knri*4+2];
     Double_t randcoinf1n=par[knri*4+1];

     //init
     Double_t N0=par[knri*4]/par[0];

     //! parent
     //! decay with 2 neutron from P2n of parent
     Double_t returnval=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf2n-randcoinfgt0n);

     //! random coinc of beta decay of parent
     returnval+=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf2n;
     //! random 1n decay of parent
     returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

     //! decay with 1 neutron from P2n of parent - randomly correlated
     returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

     for (Int_t i=0;i<npaths;i++){
         //! decay with 2 neutron from P2n of daugter
         returnval+=neueff*neueff*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(1-randcoinf2n-randcoinfgt0n);
         //! random coinc of beta decay of daugter
         returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*randcoinf2n;
         //! random 1n decay of daugter
         returnval+=neueff*pn[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

         //! decay with 1 neutron from P2n of daugter - randomly correlated
         returnval+=2*(neueff*(1-neueff))*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);
     }
     return returnval;
 }


 //! Global Bateaman function
 Double_t fitF::corefcn(Int_t ndecay,Int_t*  idecaymap,Int_t*  inneu,Double_t* production_yield, Double_t* b1n,Double_t* b2n,Double_t* lamda,Double_t N0,Double_t t) const{
     Double_t fcnret=0;

     Double_t factor1=1.;
     //only parrent decay p2n
     for (int i=0;i<ndecay-1;i++){
         Double_t corrproductionyield = production_yield[idecaymap[i+1]];

         if (flag_sum_isomer_ratio){
             for (Int_t j=0;j<nisomers;j++){
                 if (idecaymap[i+1] == groundstate[j]){
                     corrproductionyield = 1-production_yield[isomerstate[j]];
                 }
             }
         }

         if (inneu[i]==0){
             factor1=factor1 * corrproductionyield*(1-b1n[idecaymap[i]]-b2n[idecaymap[i]])*lamda[idecaymap[i]];//branching here!
         }else if (inneu[i]==1){
             factor1=factor1 * corrproductionyield*b1n[idecaymap[i]]*lamda[idecaymap[i]];
         }else{
             factor1=factor1 * corrproductionyield*b2n[idecaymap[i]]*lamda[idecaymap[i]];
         }

     }

     Double_t factor2=0;
     for (int i=0;i<ndecay;i++){
         Double_t factor2i=exp(-lamda[idecaymap[i]]*t);

         Double_t factor2ij=1;
         for (int j=0;j<ndecay;j++)
             if (j!=i) factor2ij=factor2ij*(lamda[idecaymap[j]]-lamda[idecaymap[i]]);
         factor2=factor2+factor2i/factor2ij;
     }

     fcnret=factor1*N0*factor2;
     return fcnret;
 }






