using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class WizardBoss : MonoBehaviour
{

    private enum State
    {
        Idle,
        Spawn,
        Attack,
        Hit,
        Dead
    };



    #region Variables
    public bool isDebug = false;
    public float maxHealth;
    private float health;
    private State currentState;

    private Image healthBar;
    private bool isDead;


    [Header("Direction")]
    public int damageDirection;

    [Header("Movement")]
    public int facingDirection;
    public float movementSpeed;

    [Header("EnemySpawn")]
    public Transform spawnPoint;
    private int currentMiniBoss;
    public GameObject[] miniBoss;

    public float spawnRate;
    private float spawnCD;
    private bool isMiniBossPresent;

    [Header("Attack")]
    public bool isInAttackMode = false;
    public float attackRadius;
    public Transform attackPoint;
    public float attackDamage;
    public float attackRate;
    private float attackCD;
    private Transform player;
    private bool isAttacking = false;

    [Header("Dash")]

    public GameObject enemyAfterImagePool;
    public Transform dashCheck;
    public float dashCheckRadius;
    private bool
        isDashing = false,
        canDash = true;

    [SerializeField]
    public float
        dashTime,
        dashSpeed,
        distanceBetweenImages,
        dashCoolDown;
    private float
        dashTimeLeft,
        lastImageXPos,
        lastDash = -100;

    [Header("KnockBack")]
    [SerializeField]
    private float knockBackDuration;
    [SerializeField]
    private Vector2 knockBackSpeed;
    private float
          knockBackStartTime;
    protected bool knockBegan = false;
    protected bool knockback;
    public Vector2 knockBackSpeedAttack1;
    public Vector2 knockBackSpeedAttack2;
    public Vector2 knockBackSpeedAttack3;

    [Header("Touch Hit box")]
    public Transform
      touchDamageCheck;
    private float
        lastTouchDamageTime;
    public float
        touchDamageWidth,
        touchDamageHeight;
    [SerializeField]
    private float
        touchDamageCoolDown,
        touchDamage;
    public Vector2 touchDamageBotLeft;
    public Vector2 touchDamageTopRight;

    private float[] attackDetails;

    [Header("Layer Mask")]
    [SerializeField]
    private LayerMask whatIsPlayer;
    [SerializeField]
    private LayerMask whatIsPlayerDamageable;

    private Animator anim;
    private Rigidbody2D rb;
    #endregion

    

    // Start is called before the first frame update
    void Start()
    {
        anim = GetComponent<Animator>();
        rb = GetComponent<Rigidbody2D>();
        player = GameObject.FindGameObjectWithTag("Player").transform;
        healthBar = GameObject.Find("HealthBarFill").GetComponent<Image>();
        healthBar.fillAmount = 1;
        isDead = false;
        spawnCD = spawnRate;
        currentState = State.Idle;
        currentMiniBoss = 0;
        isMiniBossPresent = false;

        health = maxHealth;
        attackDetails = new float[2];
        Instantiate(enemyAfterImagePool);
    }

    // Update is called once per frame
    void Update()
    {

        switch (currentState)
        {

            case State.Idle:
                UpdateIdleState();
                break;
            case State.Attack:
                UpdateAttackState();
                break;

        }
        
        CheckKnockBack();
        CheckDash();
        CheckTouchDamage();
    }

    #region Idle State
    private void EnterIdleState()
    {        
        spawnCD = spawnRate;
    }
    private void UpdateIdleState()
    {
        if (!isDebug)
        {
            if (spawnCD > 0 && !isMiniBossPresent)
            {
                spawnCD -= Time.deltaTime;
            }

            if (spawnCD <= 0 && currentMiniBoss < 3)
            {
                SwitchState(State.Spawn);
            }
        }

        SwitchState(State.Attack);
    }
    private void ExitIdleState()
    {

    }
    #endregion


    #region Spawn State
    private void EnterSpawnState()
    {
        anim.SetTrigger("Spawn");
        
    }
    public void SpawnEnemy()
    {
        GameObject miniBossClone = Instantiate(miniBoss[currentMiniBoss], spawnPoint.transform.position, miniBoss[currentMiniBoss].transform.rotation);
        isMiniBossPresent = true;
        SwitchState(State.Idle);
    }
    private void ExitSpawnState()
    {

    }
    #endregion

    #region Attack State
    private void EnterAttackState()
    {
        FindObjectOfType<MeleeAttackController>().BossReady();
        attackCD = attackRate;
        GameObject.Find("IceBoss").GetComponent<Renderer>().enabled = true;
        GameObject.Find("IceBoss").GetComponent<Collider2D>().enabled = true;
        GameObject.Find("MetalBoss").GetComponent<Renderer>().enabled = true;
        GameObject.Find("MetalBoss").GetComponent<Collider2D>().enabled = true;
        GameObject p = GameObject.Find("BossPad");
        p.GetComponent<Transform>().GetChild(0).gameObject.SetActive(true);
        p.GetComponent<Transform>().GetChild(1).gameObject.SetActive(true);

    }

    private void UpdateAttackState()
    {
        // Track player in order to chase if far away

        if (!isAttacking && !isDashing)
        {
            if (facingDirection == 1 && player.transform.position.x < transform.position.x)
            {
                Flip();
            }
            else if (facingDirection == -1 && player.transform.position.x > transform.position.x)
            {
                Flip();
            }


            rb.velocity = new Vector2(facingDirection * movementSpeed, rb.velocity.y);
        }

        bool playerDetected = Physics2D.OverlapCircle(attackPoint.position, attackRadius, whatIsPlayer); 

        if (attackCD > 0)
        {
            attackCD -= Time.deltaTime;
        }

        if(playerDetected)
        {
            rb.velocity = Vector2.zero;

            if (attackCD <= 0)
            {
                isAttacking = true;                
                anim.SetTrigger("Attack");
                attackCD = attackRate;
            }
        }

        bool playerDetectedBehind = Physics2D.OverlapCircle(dashCheck.position, dashCheckRadius, whatIsPlayer);

        if (playerDetectedBehind && !isAttacking && !knockback)
        {
            if (Time.time >= lastDash + dashCoolDown)
            {

                AttemptToDash();
            }

        }
        anim.SetFloat("movementSpeed", rb.velocity.magnitude);

    }
    private void CheckAttackCircle()
    {
        Collider2D detectedPlayer = Physics2D.OverlapCircle(attackPoint.position, attackRadius, whatIsPlayerDamageable);


        if (detectedPlayer)
        {
            attackDetails[0] = attackDamage;
            attackDetails[1] = transform.position.x;

            //Instantiate(hitParticle, attackPoint.transform.position, Quaternion.Euler(0, 0, Random.Range(0, 360)));
            //collider.transform.SendMessage("SetKnockBackValues", (attackComboCount - 1) % 3);

            detectedPlayer.transform.SendMessage("Damage", attackDetails);

            //FindObjectOfType<AudioManager>().Play("AttackHit");
        }

    }

    private void AttackCompleted()
    {
        isAttacking = false;
    }

    public void AttemptToDash()
    {
        isDashing = true;
        dashTimeLeft = dashTime;
        lastDash = Time.time;
        EnemyAfterImagePool.Instance.GetFromPool();
        lastImageXPos = transform.position.x;
    }

    private void CheckDash()
    {
        if (isDashing)
        {
            if (dashTimeLeft > 0)
            {
                this.gameObject.layer = 29;
                rb.velocity = new Vector2(facingDirection, 0) * dashSpeed;
                dashTimeLeft -= Time.unscaledDeltaTime;

                if (Mathf.Abs(transform.position.x - lastImageXPos) > distanceBetweenImages)
                {
                    EnemyAfterImagePool.Instance.GetFromPool();
                    lastImageXPos = transform.position.x;
                }
            }

            
            if (dashTimeLeft <= 0)
            {
                this.gameObject.layer = 11;
                isDashing = false;
                rb.velocity = Vector2.zero;
            }
        }
    }
    private void CheckTouchDamage()
    {
        if (Time.time >= lastTouchDamageTime + touchDamageCoolDown && !isDead)
        {
            touchDamageBotLeft.Set(touchDamageCheck.position.x - touchDamageWidth / 2, touchDamageCheck.position.y - touchDamageHeight / 2);
            touchDamageTopRight.Set(touchDamageCheck.position.x + touchDamageWidth / 2, touchDamageCheck.position.y + touchDamageHeight / 2);

            Collider2D hit = Physics2D.OverlapArea(touchDamageBotLeft, touchDamageTopRight, whatIsPlayerDamageable);

            if (hit != null)
            {
                hit.GetComponent<Rigidbody2D>().velocity = Vector2.zero;
                rb.velocity = Vector2.zero;
                lastTouchDamageTime = Time.time;
                attackDetails[0] = touchDamage / 2;
                attackDetails[1] = transform.position.x;
                hit.SendMessage("Damage", attackDetails);

                if (this.GetComponent<BossMage>())
                {
                    GetComponent<BossMage>().DisableChase();
                }
            }
        }
    }

    private void ExitAttackState()
    {

    }
    #endregion


    #region Hurt
    public void Damage(float[] attackDetails)
    {
        int direction;
        int visualDamage = (int)attackDetails[0];

        // currentHealth -= attackDetails[0];
        if (attackDetails[1] > transform.position.x)
        {
            damageDirection = -1;
        }
        else
        {
            damageDirection = 1;
        }

        //Instantiate(hitParticle, transform.position, Quaternion.Euler(0, 0, Random.Range(0, 360)));

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

        float healthPercentage = health / maxHealth;
        healthPercentage = Mathf.Clamp(healthPercentage, 0, 1);
        healthBar.fillAmount = healthPercentage;
        lastDash = Time.time;
        if (!isAttacking)
            anim.SetTrigger("Hit");

        if (health <= 0)
        {
            isDead = true;
            GetComponent<Collider2D>().enabled = false;
            rb.isKinematic = true;
            if (anim == null) return;
            anim.updateMode = AnimatorUpdateMode.UnscaledTime;
            anim.SetTrigger("Dead");
        }


        //SwitchState(State.Knockback);

    }
    private void onDeath()
    {
        
        Destroy(gameObject);
    }

    public void KnockBack(int direction)
    {
        knockback = true;
        knockBackStartTime = Time.time;
        rb.velocity = new Vector2(knockBackSpeed.x * direction, knockBackSpeed.y);
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
    #endregion



    private void SwitchState(State newState)
    {
        switch(currentState)
        {
            case State.Idle:
                ExitIdleState();
                break;
            case State.Spawn:
                ExitSpawnState();
                break;
            case State.Attack:
                ExitAttackState();
                break;
        }

        currentState = newState;
        switch (newState)
        {
            
            case State.Idle:                
                EnterIdleState();
                break;
            case State.Spawn:
                EnterSpawnState();
                break;
            case State.Attack:
                EnterAttackState();
                break;
        }
    }

    public void Flip()
    {
        facingDirection *= -1;
        transform.Rotate(new Vector3(0, 180.0f, 0));
    }

    public void miniBossDead()
    {
        currentMiniBoss++;
        health -= 200;

        float healthPercentage = health / maxHealth;
        healthPercentage = Mathf.Clamp(healthPercentage, 0, 1);
        healthBar.fillAmount = healthPercentage;
        isMiniBossPresent = false;

        if(currentMiniBoss >= 3)
        {
           // transform.position = GameObject.Find("TeleportPoint").transform.position;
            StartCoroutine(BeginSwitchAttackState());
        }
    }

    private IEnumerator BeginSwitchAttackState()
    {
        anim.SetTrigger("Jump");
        rb.AddForce(new Vector2(2, 1) * facingDirection * 100);
        yield return new WaitForSeconds(0.5f);
        SwitchState(State.Attack);
    }

    private void OnDrawGizmos()
    {

        Gizmos.DrawWireSphere(attackPoint.position, attackRadius);
        Gizmos.DrawWireSphere(dashCheck.position, dashCheckRadius);

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
