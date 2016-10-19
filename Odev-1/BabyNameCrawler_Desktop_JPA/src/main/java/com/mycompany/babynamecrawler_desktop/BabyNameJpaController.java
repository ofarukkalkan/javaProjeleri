/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.babynamecrawler_desktop;

import com.mycompany.babynamecrawler_desktop.exceptions.NonexistentEntityException;
import java.io.Serializable;
import java.util.List;
import javax.persistence.EntityManager;
import javax.persistence.EntityManagerFactory;
import javax.persistence.Query;
import javax.persistence.EntityNotFoundException;
import javax.persistence.criteria.CriteriaQuery;
import javax.persistence.criteria.Root;

/**
 *
 * @author gldev
 */
public class BabyNameJpaController implements Serializable {

    public BabyNameJpaController(EntityManagerFactory emf) {
        this.emf = emf;
    }
    private EntityManagerFactory emf = null;

    public EntityManager getEntityManager() {
        return emf.createEntityManager();
    }

    public void create(BabyName babyName) {
        EntityManager em = null;
        try {
            em = getEntityManager();
            em.getTransaction().begin();
            em.persist(babyName);
            em.getTransaction().commit();
        } finally {
            if (em != null) {
                em.close();
            }
        }
    }

    public void edit(BabyName babyName) throws NonexistentEntityException, Exception {
        EntityManager em = null;
        try {
            em = getEntityManager();
            em.getTransaction().begin();
            babyName = em.merge(babyName);
            em.getTransaction().commit();
        } catch (Exception ex) {
            String msg = ex.getLocalizedMessage();
            if (msg == null || msg.length() == 0) {
                Long id = babyName.getId();
                if (findBabyName(id) == null) {
                    throw new NonexistentEntityException("The babyName with id " + id + " no longer exists.");
                }
            }
            throw ex;
        } finally {
            if (em != null) {
                em.close();
            }
        }
    }

    public void destroy(Long id) throws NonexistentEntityException {
        EntityManager em = null;
        try {
            em = getEntityManager();
            em.getTransaction().begin();
            BabyName babyName;
            try {
                babyName = em.getReference(BabyName.class, id);
                babyName.getId();
            } catch (EntityNotFoundException enfe) {
                throw new NonexistentEntityException("The babyName with id " + id + " no longer exists.", enfe);
            }
            em.remove(babyName);
            em.getTransaction().commit();
        } finally {
            if (em != null) {
                em.close();
            }
        }
    }

    public List<BabyName> findBabyNameEntities() {
        return findBabyNameEntities(true, -1, -1);
    }

    public List<BabyName> findBabyNameEntities(int maxResults, int firstResult) {
        return findBabyNameEntities(false, maxResults, firstResult);
    }

    private List<BabyName> findBabyNameEntities(boolean all, int maxResults, int firstResult) {
        EntityManager em = getEntityManager();
        try {
            CriteriaQuery cq = em.getCriteriaBuilder().createQuery();
            cq.select(cq.from(BabyName.class));
            Query q = em.createQuery(cq);
            if (!all) {
                q.setMaxResults(maxResults);
                q.setFirstResult(firstResult);
            }
            return q.getResultList();
        } finally {
            em.close();
        }
    }

    public BabyName findBabyName(Long id) {
        EntityManager em = getEntityManager();
        try {
            return em.find(BabyName.class, id);
        } finally {
            em.close();
        }
    }

    public int getBabyNameCount() {
        EntityManager em = getEntityManager();
        try {
            CriteriaQuery cq = em.getCriteriaBuilder().createQuery();
            Root<BabyName> rt = cq.from(BabyName.class);
            cq.select(em.getCriteriaBuilder().count(rt));
            Query q = em.createQuery(cq);
            return ((Long) q.getSingleResult()).intValue();
        } finally {
            em.close();
        }
    }
    
}
