using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Pathfinding;

public class FlyingEnemy : Enemy
{
    float originalY;
    public float floatStrength = 1; // You can change this in the Unity Editor to change the range of y positions that are possible.      
    public AIPath aiPath; //AI Path for Pathfinder
    public Transform attackPoint; //Area to begin attacking
    public float attackRadius; //Size of Attack Radius
    public float attackRate;
    public float attackCD; //Cooldown for attack
    public float agroRange; //Agro Range for enemy
    private Transform target; //Target of enemy, typically the player
    void Start()
    {
        target = GameObject.FindGameObjectWithTag("Player").transform;
    }

    // Update is called once per frame
    void Update()
    {
        //Utilize circle to detect attack range
        base.CheckKnockBack();
        base.CheckTouchDamage();

        Collider2D detectedObjects = Physics2D.OverlapCircle(attackPoint.position, attackRadius);

        //Code for Agro 
        float distToPlayer = Vector2.Distance(transform.position, target.position);
        //If player is in agro range, begin chase. 
        if (knockback == false)
        {
            if (distToPlayer < agroRange && isDead == false)
            {
                restartAI();
            }
            //End chase when player is outside of range
            else
            {
                //stop movement to attack
                aiPath.canMove = false;
                aiPath.canSearch = false;
            }
        }

        if (Physics2D.OverlapCircle(attackPoint.position, attackRadius, whatIsPlayer) && (attackCD <= 0))
        {
            attackCD = attackRate;
            animator.SetBool("isAttacking", true);
            //stop movement to attack
            aiPath.canMove = false;
            aiPath.canSearch = false;
        }
        //Code to flip sprite according to direction of pathfinder AI
        if (aiPath.desiredVelocity.x >= 0.01f)
        {
            transform.localScale = new Vector2(-1, 1);
        }
        else if(aiPath.desiredVelocity.x <= 0.01f)
        {
            transform.localScale = new Vector2(1, 1);
        }
        if(attackCD > 0)
        {
            attackCD -= Time.deltaTime;
        }
    }

    public void restartAI()
    {
        //Restart AI movement when not attacking
        animator.SetBool("isAttacking", false);
        aiPath.canMove = true;
        aiPath.canSearch = true;
    }
}
