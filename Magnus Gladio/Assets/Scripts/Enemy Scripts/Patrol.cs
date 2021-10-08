using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine;

public class Patrol : Enemy
{
    private bool movingR; //Is sprite moving to the right
    [Space(10)]
    [Header("Patrol vars")]
    public Transform groundDetection; //Detects if sprite is making contact w/ ground
    public Transform wallDetection; //Detects when sprite has reached a wall
    public Transform attackPoint;
    public float attackRadius, attackRate;
    private float attackCD;
    private bool isWalking;
    [SerializeField]
    private LayerMask whatIsGround, whatIsWall;
    private bool playerInRange = false;
    private int facingDirection = 1;

    private GameManager manager;
    public Sound walkingAudio;


    private void Start()
    {
        isWalking = false;
        movingR = true;
        manager = GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>();
        walkingAudio.source = enemyAudioSource;
        //moveSpeed = 2f;
    }
    private void Update()
    {
        if (manager.isGameStart)
        {
            base.CheckKnockBack();
            base.CheckTouchDamage();
            base.CheckStun();

            if (!knockBegan && !playerInRange && !isStunned)
            {
                isWalking = true;
                PatrolMove();
            }
            else
            {
                isWalking = false;

                
            }
            animator.SetBool("isWalking", isWalking);

            Collider2D collider = Physics2D.OverlapCircle(attackPoint.position, attackRadius, whatIsPlayer);
            if (collider != null && collider.CompareTag("Player"))
            {
                if (attackCD <= 0)
                {
                    attackCD = attackRate;
                    playerInRange = true;
                    animator.SetBool("isAttacking", true);
                }
            }
            else
            {
                playerInRange = false;
            }

            if (attackCD > 0)
            {
                attackCD -= Time.deltaTime;
            }
        }
    }

    public void PatrolMove()
    {
        rb.velocity = new Vector2(moveSpeed * facingDirection, rb.velocity.y);

        //Play Walking Audio
        /*if (walkingAudio.source != null)
        {
            if (!walkingAudio.source.isPlaying)
            {
                enemyAudioSource.clip = walkingAudio.clip;
                enemyAudioSource.Play();
            }
        }*/

        //transform.Translate(Vector2.right * moveSpeed * Time.deltaTime); //Begin by moving character to the right

        //isGrounded = Physics2D.OverlapCircle(groundCheck.position, groundCheckRadius, whatIsGround);

        RaycastHit2D ground = Physics2D.Raycast(groundDetection.position, Vector2.down, 1.5f, whatIsGround); //Creating downward raycast in front of enemy to detect end of platform
        RaycastHit2D wall = Physics2D.Raycast(wallDetection.position, transform.right, 2f, whatIsWall); //Forward raycast to detect walls and other obstacles

        //If raycast is not detecting ground anymore, flip direction of movement
        if (ground.collider == false || wall.collider == true  )
        {            
            Flip();                                     
        }


    }

    public void Flip()
    {
        movingR = !movingR;
        transform.Rotate(0, 180, 0f, 0);
        facingDirection *= -1;
    }
    private void CheckAttackCircle()
    {
        Collider2D collider = Physics2D.OverlapCircle(attackPoint.position, attackRadius, whatIsPlayer); 
       

        if (collider != null && collider.CompareTag("Player"))
        {           
            float[] attackDetails = new float[2];

            attackDetails[0] = damage/2;
            attackDetails[1] = transform.position.x;

            collider.transform.SendMessage("Damage", attackDetails);

        }
    }
    

    private void OnDrawGizmos()
    {
        base.OnDrawGizmos();
        Gizmos.color = Color.white;
        Gizmos.DrawLine(groundDetection.position, new Vector3(groundDetection.position.x, groundDetection.position.y - 1.5f, groundDetection.position.z));
        Gizmos.DrawLine(wallDetection.position, new Vector3(wallDetection.position.x + 2, wallDetection.position.y, wallDetection.position.z));
        Gizmos.DrawWireSphere(attackPoint.position, attackRadius);
    }
}
