using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Enemy : MonoBehaviour
{
    public bool isMini; //Dictates to code whether this is full boss, or mini-boss version spawned by final boss
    public float health; //Health variable for enemy
    public float damage; //Amount of damage enemy does;
    public float moveSpeed; //Movement speed
    protected bool isDead; //Is enemy dead
    
    protected bool isStunned;
    private float stunCD;
    public float stunDuration;

    public float coinCount;
    public GameObject coinPrefab;

    protected Collider2D col; //Collider for sprite

    [SerializeField]
    protected Animator animator;

    [SerializeField]
    public Transform
      touchDamageCheck, 
        stunText;      

    [SerializeField]
    public LayerMask whatIsPlayer;

    public float[] attackDetails = new float[2];

    [Space(10)]
    [Header("KnockBack")]
    [SerializeField]
    private float knockBackDuration;
    [SerializeField]
    private Vector2 knockBackSpeed;
    private float
          knockBackStartTime;
    protected bool knockBegan = false;
    protected bool knockback;


    [Space(10)]
    [Header("Touch Damage")]
    private float
        lastTouchDamageTime;
    [SerializeField]
    public Vector2
        touchDamageBotLeft,
        touchDamageTopRight;
    public float
        touchDamageWidth,
        touchDamageHeight;   
    [SerializeField] 
    private float 
        touchDamageCoolDown,
        touchDamage;


    [Space(10)]
    [Header("Targeted")]
    public GameObject enemyTargeted;
    private GameObject enemyTargetedCone;

    protected Rigidbody2D rb;

    [Space(10)]
    [Header("Audio")]
    public AudioSource enemyAudioSource;


    protected bool isBoss = false;

    protected Image healthBar;
    protected float maxHealth;

    // Start is called before the first frame update
    private void Awake()
    {
        isDead = false;
        col = GetComponent<Collider2D>();
        animator = GetComponent<Animator>();
        if(animator != null)
            animator.SetBool("isAttacking", false);
        rb = GetComponent<Rigidbody2D>();
       
    }

    public void KnockBack(int direction)
    {
        knockback = true;
        knockBackStartTime = Time.time;
        rb.velocity = new Vector2(knockBackSpeed.x * direction, knockBackSpeed.y);
    }

    protected virtual void CheckTouchDamage()
    {
        if (Time.time >= lastTouchDamageTime + touchDamageCoolDown && !isDead)
        {
            touchDamageBotLeft.Set(touchDamageCheck.position.x - touchDamageWidth / 2, touchDamageCheck.position.y - touchDamageHeight / 2);
            touchDamageTopRight.Set(touchDamageCheck.position.x + touchDamageWidth / 2, touchDamageCheck.position.y + touchDamageHeight / 2);

            Collider2D hit = Physics2D.OverlapArea(touchDamageBotLeft, touchDamageTopRight, whatIsPlayer);

            if (hit != null)
            {
                hit.GetComponent<Rigidbody2D>().velocity = Vector2.zero;
                rb.velocity = Vector2.zero;
                lastTouchDamageTime = Time.time;
                attackDetails[0] = touchDamage/2;
                attackDetails[1] = transform.position.x;
                hit.SendMessage("Damage", attackDetails);

                if(this.GetComponent<BossMage>())
                {
                    GetComponent<BossMage>().DisableChase();
                }
            }
        }
    }
    protected virtual void CheckKnockBack()
    {
        if (Time.time >= knockBackStartTime + knockBackDuration && knockback)
        {
            knockback = false;
            knockBegan = false;
            rb.velocity = new Vector2(0.0f, rb.velocity.y);
        }

        if (lastTouchDamageTime >= 0)
        {
            lastTouchDamageTime -= Time.deltaTime;
        }
    }

    protected virtual void CheckStun()
    {
        if (isStunned)
        {
            if (stunCD > 0)
            {
                stunCD -= Time.deltaTime;
            }
            if (stunCD <= 0)
            {
                stunCD = 0;
                isStunned = false;                
            }
        }
    }

    public void finishAttack()
    {
        if(animator != null)
            animator.SetBool("isAttacking", false);
    }
    public void Damage(float[] attackDetails)
    {

        int direction;
        int visualDamage = (int)attackDetails[0];

        //controller.StopAllBasicMovement();
        GetComponent<Rigidbody2D>().velocity = Vector2.zero;

        if (attackDetails[1] < transform.position.x)
        {
            direction = 1;
        }
        else
        {
            direction = -1;
        }

        knockBegan = true;

        KnockBack(direction);      

        health -= visualDamage;

        if (isBoss && !isMini)
        {
            float healthPercentage = health / maxHealth;
            healthPercentage = Mathf.Clamp(healthPercentage, 0, 1);
            healthBar.fillAmount = healthPercentage;
        }

        if (health <= 0)
        {

            moveSpeed = 0;
            isDead = true;

            GetComponent<BoxCollider2D>().enabled = false;
            GetComponent<CircleCollider2D>().enabled = false;
            this.gameObject.layer = 8;
            GetComponent<Rigidbody2D>().isKinematic = true;

            DropCoins();

            if (animator == null) return;

            animator.updateMode = AnimatorUpdateMode.UnscaledTime;
            animator.SetBool("isDead", true);
           
        }


        
    }

    public void Stun()
    {
        GetComponent<Rigidbody2D>().velocity = Vector2.zero;
        isStunned = true;
        stunCD = stunDuration;
        Vector3 pos = new Vector3(0, GetComponent<BoxCollider2D>().size.y + 1, 0f);
        Instantiate(stunText,transform.position + pos, Quaternion.identity);
        Destroy(enemyTargetedCone);
    }

    public void Targeted()
    {
        enemyTargetedCone = Instantiate(enemyTargeted, transform);
    }
    public void UnTarget()
    {
        if(enemyTargetedCone != null)
        {
            Destroy(enemyTargetedCone);
        }
    }
   
    private void DropCoins()
    {
        for(int i = 0; i < coinCount; i++)
        {
            GameObject coinClone = Instantiate(coinPrefab, transform.position, Quaternion.identity);
            Vector2 randomDirection = new Vector3(Random.value, Random.value);
            coinClone.GetComponent<Rigidbody2D>().AddForce(randomDirection.normalized * 200);
        }

    }
    private void onDeath()
    {
        //GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>().InstantiateHealthPotion(transform.localPosition);

        if(isBoss && !isMini)
        {
            GameObject.FindGameObjectWithTag("Teleporter").GetComponent<Teleporter>().EnableTeleporter();
        }
        
        if(isMini)
        {
            GameObject.FindObjectOfType<WizardBoss>().miniBossDead();
        }

        Destroy(gameObject);
    }

    public void OnDrawGizmos()
    {
        Gizmos.color = Color.black;
        Vector2 botLeft = new Vector2(touchDamageCheck.position.x - touchDamageWidth / 2, touchDamageCheck.position.y - touchDamageHeight / 2);
        Vector2 botRight = new Vector2(touchDamageCheck.position.x + touchDamageWidth / 2, touchDamageCheck.position.y - touchDamageHeight / 2);
        Vector2 topLeft = new Vector2(touchDamageCheck.position.x - touchDamageWidth / 2, touchDamageCheck.position.y + touchDamageHeight / 2);
        Vector2 topRight = new Vector2(touchDamageCheck.position.x + touchDamageWidth / 2, touchDamageCheck.position.y + touchDamageHeight / 2);
        Gizmos.DrawLine(botLeft, botRight);
        Gizmos.DrawLine(botRight, topRight);
        Gizmos.DrawLine(topRight, topLeft);
        Gizmos.DrawLine(topLeft, botLeft);
    }
}
