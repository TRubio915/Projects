using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine;

public class ChaseEnemy : Enemy
{
    public Transform target; //Target of enemy, typically the player
    public float agroRange; //Agro Range for enemy

    private Rigidbody2D rb2d; //Sprites rigidbody

    public Transform groundDetection; //Detects if sprite is making contact w/ ground
    public Transform wallDetection; //Detects when sprite has reached a wall

    [SerializeField]
    private LayerMask whatIsGround, whatIsWall; //Layer masks so collision knows when it is in contact w/ ground or wall

    [SerializeField]
    private float groundCheckDistance, wallCheckDistance;

    private float startYPos;

    public Sound soundScreech;


    void Start()
    {
        rb2d = GetComponent<Rigidbody2D>();
        target = GameObject.FindGameObjectWithTag("Player").transform;
        startYPos = transform.position.y;
        soundScreech.source = enemyAudioSource;
    }

    // Update is called once per frame
    void Update()
    {
        //Check distance to player
        CheckKnockBack();
        CheckTouchDamage();

        float distToPlayer = Vector2.Distance(transform.position, target.position);
        //If player is in agro range, begin chase. 
        if (knockback == false)
        {
            if (distToPlayer < agroRange && isDead == false)
            {
                ChasePlayer();
            }
            //End chase when player is outside of range
            else
            {
                StopChase();
            }
        }
    }

    protected override void CheckKnockBack()
    {
        base.CheckKnockBack();
        Vector3 newPos = new Vector3(transform.position.x, startYPos, transform.position.z);
        transform.position = Vector3.Lerp(transform.position, newPos, 0.5f);
    }

    private void ChasePlayer()
    {
        animator.SetBool("isChasing", true); //Start chase animation
        RaycastHit2D ground = Physics2D.Raycast(groundDetection.position, Vector2.down, groundCheckDistance, whatIsGround); //Creating downward raycast in front of enemy to detect end of platform
        RaycastHit2D wall = Physics2D.Raycast(wallDetection.position, Vector2.left, wallCheckDistance, whatIsWall); //Forward raycast to detect walls and other obstacles

        //Play Ghost Chase Audio
        
        if (soundScreech.source != null)
        {
            if (soundScreech.source.isPlaying == false)
            {
                enemyAudioSource.clip = soundScreech.clip;
                enemyAudioSource.Play();
            }
        }
        //If raycast is not detecting ground anymore or hits wall stop chase
        if (ground.collider == false || wall.collider == true)
        {
            StopChase();
        }
        //Decide movement direction
        if (transform.position.x < target.position.x)
        {
            //enemy left of target, move right
            rb2d.velocity = new Vector2(moveSpeed, 0);
            transform.localScale = new Vector2(-1, 1);
        }
        else
        {
            //enemy right of target, move left
            rb2d.velocity = new Vector2(-moveSpeed, 0);
            transform.localScale = new Vector2(1, 1);
        }
    }

    private void StopChase()
    {
        //Stop Sound
        FindObjectOfType<AudioManager>().Stop("Ghost");
        //Stop movement
        rb2d.velocity = Vector2.zero;
        //Stop chase animation
        animator.SetBool("isChasing", false);
    }

    private void OnDrawGizmos()
    {
        base.OnDrawGizmos();
        //Draw lines so detection is more clear
        Gizmos.DrawLine(groundDetection.position, new Vector3(groundDetection.position.x, groundDetection.position.y - groundCheckDistance, groundDetection.position.z));
        Gizmos.DrawLine(wallDetection.position, new Vector3(wallDetection.position.x - wallCheckDistance, wallDetection.position.y, wallDetection.position.z));
    }
}
