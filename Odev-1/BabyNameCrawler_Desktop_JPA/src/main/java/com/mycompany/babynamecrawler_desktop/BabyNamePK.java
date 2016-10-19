/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.babynamecrawler_desktop;

import java.io.Serializable;
import javax.persistence.Basic;
import javax.persistence.Column;
import javax.persistence.Embeddable;
import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

/**
 *
 * @author gldev
 */
@Embeddable
public class BabyNamePK implements Serializable {

    @Basic(optional = false)
    @NotNull
    @Column(name = "YEARR")
    private int yearr;
    @Basic(optional = false)
    @NotNull
    @Size(min = 1, max = 50)
    @Column(name = "NAME")
    private String name;
    @Basic(optional = false)
    @NotNull
    @Column(name = "GENDER")
    private Character gender;

    public BabyNamePK() {
    }

    public BabyNamePK(int yearr, String name, Character gender) {
        this.yearr = yearr;
        this.name = name;
        this.gender = gender;
    }

    public int getYearr() {
        return yearr;
    }

    public void setYearr(int yearr) {
        this.yearr = yearr;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public Character getGender() {
        return gender;
    }

    public void setGender(Character gender) {
        this.gender = gender;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (int) yearr;
        hash += (name != null ? name.hashCode() : 0);
        hash += (gender != null ? gender.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BabyNamePK)) {
            return false;
        }
        BabyNamePK other = (BabyNamePK) object;
        if (this.yearr != other.yearr) {
            return false;
        }
        if ((this.name == null && other.name != null) || (this.name != null && !this.name.equals(other.name))) {
            return false;
        }
        if ((this.gender == null && other.gender != null) || (this.gender != null && !this.gender.equals(other.gender))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "com.mycompany.babynamecrawler_desktop.BabyNamePK[ yearr=" + yearr + ", name=" + name + ", gender=" + gender + " ]";
    }
    
}
