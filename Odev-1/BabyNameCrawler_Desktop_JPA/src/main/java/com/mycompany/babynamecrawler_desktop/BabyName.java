/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.babynamecrawler_desktop;

import java.io.Serializable;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;

/**
 *
 * @author gldev
 */
@Entity
public class BabyName implements Serializable {

    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy = GenerationType.AUTO)
    private Long id;
    
    /*benim eklediklerim*/
    private int yearr;
    private String name;
    private char gender;
    private int count;

    public void setYearr(int yearr) {
        this.yearr = yearr;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setGender(char gender) {
        this.gender = gender;
    }

    public void setCount(int count) {
        this.count = count;
    }

    public int getYearr() {
        return yearr;
    }

    public String getName() {
        return name;
    }

    public char getGender() {
        return gender;
    }

    public int getCount() {
        return count;
    }

    public Long getId() {
        return id;
    }

    public void setId(Long id) {
        this.id = id;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (id != null ? id.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BabyName)) {
            return false;
        }
        BabyName other = (BabyName) object;
        if ((this.id == null && other.id != null) || (this.id != null && !this.id.equals(other.id))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "com.mycompany.babynamecrawler_desktop.BabyName[ id=" + id + " ]";
    }
    
}
