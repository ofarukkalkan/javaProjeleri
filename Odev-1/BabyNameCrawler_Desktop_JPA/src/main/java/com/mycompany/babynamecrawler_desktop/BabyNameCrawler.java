/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.babynamecrawler_desktop;


import java.net.*;
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.persistence.EntityManagerFactory;
import javax.persistence.Persistence;


/**
 *
 * @author gldev
 */
public class BabyNameCrawler {

    public static void main(String[] args) {
        String PERSISTENCE_UNIT_NAME = "BabyName";
        EntityManagerFactory factory =factory = Persistence.createEntityManagerFactory(PERSISTENCE_UNIT_NAME);

        BabyNameJpaController kontrolcu = new BabyNameJpaController(factory);
        
        BabyName name = new BabyName(); // cekilecek isimler jpa ile eslesen bu objede tutulacak
        BabyNamePK pk = new BabyNamePK();
        

        List<String> liste = new ArrayList<String>();

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
            System.out.println(mue.toString());
        } catch (IOException ioe) {
            System.out.println(ioe.toString());
        } finally {
            try {
                if (is != null) {
                    is.close();
                }
            } catch (IOException ioe) {
            }
        }
        for (int i = 0; i < 1000; i++) {

            kelimeler = liste.get(i).replaceAll(" ", "").split("\t");
            
            pk.setYearr(2002);
            pk.setName(kelimeler[1]);
            pk.setGender('M');
            name.setCount(Integer.parseInt(kelimeler[2]));
            name.setBabyNamePK(pk);
            
            try {
                kontrolcu.create(name); // erkek ismini ekle
            } catch (Exception ex) {
                Logger.getLogger(BabyNameCrawler.class.getName()).log(Level.SEVERE, null, ex);
            }
            
            
                     
            pk.setYearr(2002);
            pk.setName(kelimeler[3]);
            pk.setGender('F');
            name.setCount(Integer.parseInt(kelimeler[4]));
            name.setBabyNamePK(pk);
            
            try {
                kontrolcu.create(name); // kiz ismini ekle
            } catch (Exception ex) {
                Logger.getLogger(BabyNameCrawler.class.getName()).log(Level.SEVERE, null, ex);
            }
            
            sayac += 2;
            

        }
        System.out.println("Eklenen kayit sayisi = " + sayac);

    }
}
