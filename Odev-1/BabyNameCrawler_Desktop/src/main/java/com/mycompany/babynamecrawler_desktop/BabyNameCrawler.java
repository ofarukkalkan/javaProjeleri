package com.mycompany.babynamecrawler_desktop;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;
import java.util.ArrayList;

/**
 *
 * @author gldev
 */
public class BabyNameCrawler {

    public static void main(String[] args) {
        // baglanti icin gerekli giris bilgileri
        String host = "jdbc:derby://localhost:1527/BabyNameRanking";
        String uName = "test";
        String uPass = "123456";
        // veritabanina gonderilecek sorguyu tutacak nesne
        Statement statement = null;

        try {
            // veritabanina baglantinin kurulmasi
            Connection con = DriverManager.getConnection(host, uName, uPass);
            // baglanti kurulan veritabanindan sorgu sablonunun alinmasi
            statement = con.createStatement();
        } catch (Exception ex) {
            System.out.println(ex.toString());
        } finally {
            System.out.println("Veritabanina baglanti basarili");
        }

        ArrayList<URL> urls = new ArrayList<>(10);

        // url den veri getirmek icin gerekli nesneler
        InputStream is = null;
        BufferedReader br;
        // her bir satirin kopyalanacagi string
        String line;
        // her bir satirdan elde edilen kelime parcaciklarinin tutulmasi icin dizi
        String[] kelimeler = new String[5];
        // eklenen kayit icin sayac
        int sayac = 0;

        try {

            // 2001 den 2010 a kadar olan linklerin eklenmesi
            urls.add(new URL("http://www.cs.armstrong.edu/liang/data/babynamesranking2001.txt"));
            urls.add(new URL("http://www.cs.armstrong.edu/liang/data/babynamesranking2002.txt"));
            urls.add(new URL("http://www.cs.armstrong.edu/liang/data/babynamesranking2003.txt"));
            urls.add(new URL("http://www.cs.armstrong.edu/liang/data/babynamesranking2004.txt"));
            urls.add(new URL("http://www.cs.armstrong.edu/liang/data/babynamesranking2005.txt"));
            urls.add(new URL("http://www.cs.armstrong.edu/liang/data/babynamesranking2006.txt"));
            urls.add(new URL("http://www.cs.armstrong.edu/liang/data/babynamesranking2007.txt"));
            urls.add(new URL("http://www.cs.armstrong.edu/liang/data/babynamesranking2008.txt"));
            urls.add(new URL("http://www.cs.armstrong.edu/liang/data/babynamesranking2009.txt"));
            urls.add(new URL("http://www.cs.armstrong.edu/liang/data/babynamesranking2010.txt"));

            // dis dongu her bir yil icin calisir
            for (int i = 0; i < urls.size(); i++) {
                // satir parcalama isinde kullanilacak dinamik dizi icin
                // liste her bir yili eklemeden once sifirlanir.
                ArrayList<String> liste = new ArrayList<String>();
                System.out.println(urls.get(i)+" => okunuyor");
                is = urls.get(i).openStream();  // throws an IOException
                br = new BufferedReader(new InputStreamReader(is));
                // url den okunan her satirin listeye eklenmesi
                for (line = br.readLine(); line != null; line = br.readLine()) {
                    liste.add(line);
                }

                // bir yilda yer alan isimler icin calisir
                for (int j = 0; j < liste.size(); j++) {
                    // elde edilen satirdan kelime parcalarinin ayirt edilme islemi
                    kelimeler = liste.get(j).replaceAll(" ", "").split("\t");

                    try {
                        // kayitlarin veritabanina insert edilme islemi
                        statement.executeUpdate("INSERT INTO BabyName " + "VALUES (" + (int)(2001 + i) + ", '" + kelimeler[1] + "', 'M', " + Integer.parseInt(kelimeler[2]) + ")");
                        statement.executeUpdate("INSERT INTO BabyName " + "VALUES (" + (int)(2001 + i) + ", '" + kelimeler[3] + "', 'F', " + Integer.parseInt(kelimeler[4]) + ")");
                        sayac += 2;
                    } catch (Exception ex) {
                        System.out.println(ex.toString());
                        break;
                    }

                }
            }

        } catch (Exception ex) {
            System.out.println(ex.toString());
        } finally {

            try {
                if (is != null) {
                    is.close(); //inputStream Kapatma
                }
            } catch (Exception ex) {
                System.out.println(ex.toString());
            }
        }

        System.out.println("Eklenen kayit sayisi = " + sayac);
    }
}
