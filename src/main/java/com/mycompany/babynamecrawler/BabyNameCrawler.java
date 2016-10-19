/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.babynamecrawler;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;


/**
 *
 * @author gldev
 */
public class BabyNameCrawler {

    public static void main(String[] args) {
        String host = "jdbc:derby://localhost:1527/BabyNameRanking";
        String uName = "test";
        String uPass = "123456";
        Statement statement = null;

        try {
            Connection con = DriverManager.getConnection(host, uName, uPass);
            statement = con.createStatement();
        } catch (SQLException ex) {
            System.out.println(ex.toString());

        } finally {
            System.out.println("veritabani baglantisi basarili");

        }

        List<String> liste = new ArrayList<String>();
        List<Dugum> dugumListesi = new ArrayList<Dugum>();

        URL url;
        InputStream is = null;
        BufferedReader br;
        String line;
        String[] kelimeler = new String[1001];
        int sayac = 0;
        try {
            url = new URL("http://www.cs.armstrong.edu/liang/data/babynamesranking2002.txt");
            is = url.openStream();  // throws an IOException
            br = new BufferedReader(new InputStreamReader(is));

            for (line = br.readLine(); line != null; line = br.readLine()) {

                liste.add(line);
            }
        } catch (MalformedURLException mue) {
            mue.printStackTrace();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        } finally {
            try {
                if (is != null) {
                    is.close();
                }
            } catch (IOException ioe) {
            }
        }
        for (int i = 0; i < 1000; i++) {

            Dugum gecici = new Dugum(null, null, 0, 0);
            // kelimeler = liste.get(i).replaceAll("[0-9]", "").split("\\s+");
            kelimeler = liste.get(i).replaceAll(" ", "").split("\t");

            gecici.isimErkek = kelimeler[1];
            gecici.isimKadin = kelimeler[2];
            dugumListesi.add(gecici);
            try {
                statement.executeUpdate("INSERT INTO BabyName " + "VALUES (" + 2002 + ", '" + kelimeler[1] + "', 'M', " + Integer.parseInt(kelimeler[2]) + ")");
                statement.executeUpdate("INSERT INTO BabyName " + "VALUES (" + 2002 + ", '" + kelimeler[1] + "', 'F', " + Integer.parseInt(kelimeler[2]) + ")");

            } catch (SQLException ex) {
                System.out.println(ex.toString());
                break;
            } finally {
                sayac += 2;// her seferde kadın-erkek olarak iki kayıt ekleniyor
            }

        }
        System.out.println("Veritabanina eklenen kayit sayisi = "+sayac);
        // kelimeler=liste.get(1).split(" ");
        /*
        System.out.println(" \n1. erkek= " + dugumListesi.get(0).isimErkek + " \n1. kadim= " + dugumListesi.get(0).isimKadin);
        System.out.println(" \n2. erkek= " + dugumListesi.get(1).isimErkek + " \n2. kadin= " + dugumListesi.get(1).isimKadin);
        System.out.println(" \n3. erkek= " + dugumListesi.get(2).isimErkek + " \n3. kadin= " + dugumListesi.get(2).isimKadin);
         */
 /*for (int i = 0; i < 100; i++) {
            if(kelimeler[i]!=null)
            System.out.println(i+". kelime= "+kelimeler[i]);
        }
         */
    }
}
