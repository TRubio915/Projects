using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine;
using UnityEngine.UI;

public class IceMage : Enemy
{
    public Transform target; //Target of enemy, typically the player
    public float fireRange; //Range where enemy will perform ranged attack
    public float fireCD; //Cooldown for fire
    public float fireAttackCircleRadius;
    public float attackRate;

    //fireball variables
    public GameObject fireball; //Fireball prefab
    public float fbSpeed; //Speed of projectile
    public Transform firePoint; //Fire point of fireball

    private Rigidbody2D rb2d; //Sprites rigidbody

    private bool isAttacking = false;

    private float startYPos;

    void Start()
    {
        rb2d = GetComponent<Rigidbody2D>();
        target = GameObject.FindGameObjectWithTag("Player").transform;
        maxHealth = health;
        isBoss = false;
        startYPos = transform.position.y;
    }

    // Update is called once per frame
    void Update()
    {
        //Check distance to player
        CheckKnockBack();
        CheckTouchDamage();

        float distToPlayer = Vector2.Distance(transform.position, target.position);
        //If player is in aggro range, begin chase. 
        if (knockback == false)
        {
            //launch fireballs
            if (distToPlayer < fireRange && isDead == false)
            {
                fireAttack();
            }
        }

        //Check and adjust cooldown
        if (fireCD > 0)
        {
            fireCD -= Time.deltaTime;
        }
    }

    protected override void CheckKnockBack()
    {
        base.CheckKnockBack();
        Vector3 newPos = new Vector3(transform.position.x, startYPos, transform.position.z);
        transform.position = Vector3.Lerp(transform.position, newPos, 0.5f);
    }

    public void EndAttack()
    {
        isAttacking = false;
        animator.SetBool("isAttacking", false);
    }

    private void fireAttack()
    {
        if ((fireCD <= 0))
        {
            isAttacking = true;
            animator.SetBool("isAttacking", isAttacking);
            Instantiate(fireball, firePoint.position, firePoint.rotation);
            fireCD = attackRate / 2;
        }
    }

    private void OnDrawGizmos()
    {
        base.OnDrawGizmos();
        Gizmos.color = Color.white;
        Gizmos.DrawWireSphere(transform.position, fireRange);
    }
}
