using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine;
using UnityEngine.UI;

public class BossDemon : Enemy
{
    public Transform target; //Target of enemy, typically the player
    public float aggroRange; //Aggro Range for enemy
    public float fireRange; //Range where enemy will switch to ranged attack
    public float attackCD; //Cooldown for attack
    public float fireCD; //Cooldown for fire
    public Transform attackPoint;
    public float attackRadius; //Size of Attack Radius
    public float fireAttackCircleRadius;
    public float attackRate;
    public int faceLeft;

    //fireball variables
    public GameObject fireball; //Fireball prefab
    public float fbSpeed; //Speed of projectile
    public Transform firePoint; //Fire point of fireball

    private Rigidbody2D rb2d; //Sprites rigidbody

    public Transform groundDetection; //Detects if sprite is making contact w/ ground
    public Transform wallDetection; //Detects when sprite has reached a wall

    [SerializeField]
    private LayerMask whatIsGround, whatIsWall; //Layer masks so collision knows when it is in contact w/ ground or wall

    [SerializeField]
    private float groundCheckDistance, wallCheckDistance;

    private float startYPos;


    private bool isAttackingClose = false, isAttackingFar = false;

    void Start()
    {
        rb2d = GetComponent<Rigidbody2D>();
        target = GameObject.FindGameObjectWithTag("Player").transform;
        startYPos = transform.position.y;
        faceLeft = 1;
        maxHealth = health;
        isBoss = true;
        attackCD = attackRate;

        if (!isMini)
        {
            healthBar = GameObject.Find("HealthBarFill").GetComponent<Image>();
            healthBar.fillAmount = 1;
        }
    }
    // Update is called once per frame
    void Update()
    {
        //Check distance to player
        CheckKnockBack();
        CheckTouchDamage();

        //Collider2D detectedObjects = Physics2D.OverlapCircle(attackPoint.position, attackRadius);

        float distToPlayer = Vector2.Distance(transform.position, target.position);
        //If player is in aggro range, begin chase. 
        if (knockback == false)
        {
            //Chase and close-range attack
            if (distToPlayer < aggroRange && isDead == false && attackCD <= 0 && isMini == false)
            {
                ChasePlayer();
            }
            //Stand still and launch fireballs
            if (distToPlayer > aggroRange && distToPlayer < fireRange && isDead == false && attackCD <= 0)
            {
                fireAttack();
            }
            //Stop AI when player is outside of range
            else
            {
                StopChase();
            }
        }

        Collider2D playerColl = Physics2D.OverlapCircle(transform.position, attackRadius, whatIsPlayer);

        if ((playerColl != null && playerColl.CompareTag("Player")) && (attackCD <= 0) && !isAttackingClose)
        {
            StopChase(); //Stopchase to attack
            attackCD = attackRate;
            isAttackingClose = true;
            animator.SetBool("IsAttackingFire", isAttackingClose);
        }
        //If cooldown is active move away from player
        else if (attackCD > 0 && fireCD > 0)
        {
            transform.position = Vector2.MoveTowards(transform.position, target.position, -1 * moveSpeed * Time.deltaTime);
        }

        if ((transform.position.x < target.position.x && faceLeft == 1) || (transform.position.x > target.position.x && faceLeft == -1))
        {
            //enemy left of player, move right
            transform.Rotate(0f, 180f, 0f);
            faceLeft = faceLeft * -1;
        }

        //Check and adjust cooldowns
        if (attackCD > 0)
        {
            attackCD -= Time.deltaTime;
        }
        if (fireCD > 0)
        {
           fireCD -= Time.deltaTime;
        }
    }

    private void CheckAttackCircle()
    {
        Collider2D collider = Physics2D.OverlapCircle(attackPoint.position, fireAttackCircleRadius, whatIsPlayer);


        if (collider != null && collider.CompareTag("Player"))
        {
            float[] attackDetails = new float[2];

            attackDetails[0] = damage/2;
            attackDetails[1] = transform.position.x;

            collider.transform.SendMessage("Damage", attackDetails);

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
        isAttackingClose = false;
        animator.SetBool("IsAttackingFire", false);
        isAttackingFar = false;
        animator.SetBool("isAttacking", false);
    }

    private void ChasePlayer()
    {
        transform.position = Vector2.MoveTowards(transform.position, target.position, 1 * moveSpeed * Time.deltaTime);
        //Decide movement direction
        if ((transform.position.x < target.position.x && faceLeft == 1) || (transform.position.x > target.position.x && faceLeft == -1))
        {
            //enemy left of player, move right
            transform.Rotate(0f, 180f, 0f);
            faceLeft = faceLeft * -1;
        }
    }

    private void StopChase()
    {
        //Stop movement
        rb2d.velocity = Vector2.zero;
        //Stop chase animation
        //animator.SetBool("isChasing", false);
    }

    private void fireAttack()
    {
        if ((fireCD <= 0))
        {
            isAttackingFar = true;
            animator.SetBool("isAttacking", isAttackingFar);
            Instantiate(fireball, firePoint.position, firePoint.rotation);
            fireCD = attackRate / 2;
        }
    }

    private void OnDrawGizmos()
    {
        base.OnDrawGizmos();
        //Draw lines so detection is more clear
        Gizmos.DrawLine(groundDetection.position, new Vector3(groundDetection.position.x, groundDetection.position.y - groundCheckDistance, groundDetection.position.z));
        Gizmos.DrawLine(wallDetection.position, new Vector3(wallDetection.position.x - wallCheckDistance, wallDetection.position.y, wallDetection.position.z));

        Gizmos.DrawWireSphere(transform.position, aggroRange);

        Gizmos.color = Color.white;
        Gizmos.DrawWireSphere(transform.position, fireRange);
        Gizmos.color = Color.red;
        Gizmos.DrawWireSphere(transform.position, attackRadius);
        Gizmos.color = Color.blue;
        Gizmos.DrawWireSphere(attackPoint.position, fireAttackCircleRadius);


    }
}
