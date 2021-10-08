using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class GolemBoss : MonoBehaviour
{
    struct RandomSelection
    {
        private int minValue;
        private int maxValue;
        public float probability;

        public RandomSelection(int minValue, int maxValue, float probability)
        {
            this.minValue = minValue;
            this.maxValue = maxValue;
            this.probability = probability;
        }

        public int GetValue() { return Random.Range(minValue, maxValue + 1); }
    }

    private enum State
    {       
        Idle,
        Knockback,
        MeleeAttack,
        RangeAttack,
        LaserAttack,
        ReCharge,
        Teleport,
        Enrage,
        Dead
    }

    private State currentState;

    #region variables

    public bool isMiniBoss;
    public float health;
    public bool isDead;
    public Transform[] positions;
    private int currentPos;

    private float teleportCD;
    public float teleportRate;

    private float currentDamageTaken;

    public int currentDirection;

    private Image healthBar;
    private float maxHealth;

    [Header ("Direction")]
    public int
       facingDirection,
       damageDirection;

    [Header ("Health")]
    private float MaxHealth;
    private float currentHealth;

    [Header ("Melee Attack")]
    [SerializeField]
    private Transform attackPoint;
    private float attackCD;

    public float
        attackRate,
        meleeAttackRadius,
        meleeAttackDamage;

    private bool isMeleeAttacking;


    [Header("Range Attack")]

    [SerializeField]
    private Transform rangeAttackPoint;

    public GameObject armPrefab;
   
    private float rangeAttackCD;

    public float
        rangeAttackRate,
        rangeAttackDistance,
        rangeAttackDamage;

    private bool isRangeAttacking;

    [Header("Laser Attack")]
    [SerializeField]
    private Transform laserLaunchPoint;
    [SerializeField]
    private GameObject laserPrefab;

    private float laserAttackCD;

    public float
        laserAttackRate,
        laserAttackDistance,
        laserAttackDamage;

    [Header("Recharge")]

    private float rechargeCD;

    public float rechargeRate;

    [Header("Enrage")]
    public float enrageDuration;
    private float enrageCD;
    private float enrageTime;
    private bool isEnraged = false;
    public Transform enragePosition;

    public GameObject snowBallPrefab;

    public Transform[] snowBallSpawnPoints;

    public float snowBallInterval;


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

    [Header ("Layer Mask")]
    [SerializeField]
    private LayerMask whatIsPlayer;

    private Animator anim;
    private Rigidbody2D rb;

    #endregion

    // Start is called before the first frame update
    void Start()
    {
        rb = GetComponent<Rigidbody2D>();
        attackDetails = new float[2];
        anim = GetComponent<Animator>();

        currentDamageTaken = 0;

        laserAttackCD = 0;
        rangeAttackCD = 0;

        currentPos = 0;

        if (!isMiniBoss)
        {
            positions[0] = GameObject.Find("GolemBossPos1").transform;
            positions[1] = GameObject.Find("GolemBossPos2").transform;

            enragePosition = GameObject.Find("GolemEnragePos").transform;

            snowBallSpawnPoints[0] = GameObject.Find("SnowBallSpawnPoint1").transform; ;
            snowBallSpawnPoints[1] = GameObject.Find("SnowBallSpawnPoint2").transform; ;
            healthBar = GameObject.Find("HealthBarFill").GetComponent<Image>();
            healthBar.fillAmount = 1;
            transform.position = positions[currentPos].position;
        }

       
        teleportCD = teleportRate;

        currentDirection = -1;

        
        maxHealth = health;

      SwitchState(State.Idle);
    }

    // Update is called once per frame
    void Update()
    {
        switch (currentState)
        {
            case State.Idle:
                UpdateIdleState();
                break;

            case State.Knockback:
                //UpdateKnockBackState();
                break;

            case State.MeleeAttack:
                UpdateMeleeAttackState();
                break;
            case State.RangeAttack:
                break;

            case State.LaserAttack:
                break;

            case State.ReCharge:
                UpdateReChargeState();
                break;

            case State.Enrage:

                break;
            
            case State.Dead:
                //UpdateDeadState();
                break;

        }
        CheckTouchDamage();
    }

    private void CheckTouchDamage()
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

        if (!isMiniBoss)
        {
            float healthPercentage = health / maxHealth;
            healthPercentage = Mathf.Clamp(healthPercentage, 0, 1);
            healthBar.fillAmount = healthPercentage;
        }
        if (health <= 0)
        {
            //moveSpeed = 0;
            isDead = true;

            GetComponent<Collider2D>().enabled = false;
            rb.isKinematic = true;
            
            if (anim == null) return;

            anim.updateMode = AnimatorUpdateMode.UnscaledTime;
            anim.SetTrigger("isDead");

        }


        currentDamageTaken += visualDamage;

       

        //SwitchState(State.Knockback);

    }
    private void onDeath()
    {
        //GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>().InstantiateHealthPotion(transform.localPosition);

        if(!isMiniBoss)
            GameObject.FindGameObjectWithTag("Teleporter").GetComponent<Teleporter>().EnableTeleporter();
        else
        {
            GameObject.FindObjectOfType<WizardBoss>().miniBossDead();
        }
       
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

    #region IdleState
    private void EnterIdleState()
    {

        //Reset idle variables here
        if (transform.position.x > FindObjectOfType<PlayerController>().transform.position.x && currentDirection == 1)
        {
            currentDirection = -1;
            transform.Rotate(0, 180, 0f, 0);
        }
        else if (transform.position.x < FindObjectOfType<PlayerController>().transform.position.x && currentDirection == -1)
        {
            currentDirection = 1;
            transform.Rotate(0, 180, 0f, 0);
        }

    }
    private void UpdateIdleState()
    {
        //Make calculation for attack / defense here





        bool playerDetected = Physics2D.Raycast(attackPoint.position, transform.right, meleeAttackRadius, whatIsPlayer);

        if (playerDetected && attackCD <= 0)
        {
            SwitchState(State.MeleeAttack);
        }

        if(attackCD > 0)
        {
            attackCD -= Time.deltaTime;
        }


        bool laserRangeDetected = Physics2D.Raycast(attackPoint.position, transform.right, laserAttackDistance, whatIsPlayer);

        bool rangeAttackDetected = Physics2D.Raycast(attackPoint.position, transform.right, rangeAttackDistance, whatIsPlayer);


        if ((rangeAttackDetected && rangeAttackCD <= 0) && (laserRangeDetected && laserAttackCD <= 0) && !playerDetected)
        {

            int random = GetRandomValue(
                    new RandomSelection(0, 5, .75f),
                    new RandomSelection(6, 8, .25f));
            
            if (isMiniBoss)
                random = 5;

            if (random > 5)
            {
                SwitchState(State.LaserAttack);
            }
            else
                SwitchState(State.RangeAttack);
        }


        if (laserAttackCD > 0)
        {
            laserAttackCD -= Time.deltaTime;
        }


        if (rangeAttackCD > 0)
        {
            rangeAttackCD -= Time.deltaTime;
        }


        if(teleportCD > 0)
        {
            teleportCD -= Time.deltaTime;
        }
        if(teleportCD <= 0 && !isMiniBoss)
        {
            teleportCD = teleportRate;
            SwitchState(State.Teleport);
        }

        if (currentDamageTaken > 60 && !isMiniBoss)
        {

            currentDamageTaken = 0;

            SwitchState(State.Enrage);

        }

    }
    int GetRandomValue(params RandomSelection[] selections)
    {
        float rand = Random.value;
        float currentProb = 0;
        foreach (var selection in selections)
        {
            currentProb += selection.probability;
            if (rand <= currentProb)
                return selection.GetValue();
        }

        //will happen if the input's probabilities sums to less than 1
        //throw error here if that's appropriate
        return -1;
    }
    private void ExitIdleState()
    {
        //Set New variables here if any
    }
    #endregion

    #region MeleeAttack
    private void EnterMeleeAttackState()
    {
        //Reset idle variables here
        attackCD = attackRate;
        isMeleeAttacking = false;
    }
    private void UpdateMeleeAttackState()
    {

        if (!isMeleeAttacking)
        {
            anim.SetTrigger("isMeleeAttack");
            isMeleeAttacking = true;
        }
        
    }

    public void checkAttackCircle()
    {
        Collider2D detectedPlayer = Physics2D.OverlapCircle(attackPoint.position, meleeAttackRadius, whatIsPlayer);


        if(detectedPlayer)
        {
            attackDetails[0] = meleeAttackDamage;
            attackDetails[1] = transform.position.x;

            //Instantiate(hitParticle, attackPoint.transform.position, Quaternion.Euler(0, 0, Random.Range(0, 360)));
            //collider.transform.SendMessage("SetKnockBackValues", (attackComboCount - 1) % 3);

            detectedPlayer.transform.SendMessage("Damage", attackDetails);
            
            //FindObjectOfType<AudioManager>().Play("AttackHit");
        }
        SwitchState(State.Idle);
    }

    private void ExitMeleeAttackState()
    {
        isMeleeAttacking = false;
    }

    #endregion

    #region Laser Attack

    private void EnterLaserAttackState()
    {
        anim.SetTrigger("isArmorBuff");
    }

    public void BeginLaserAttack()
    {
        Instantiate(laserPrefab, laserLaunchPoint.position, Quaternion.identity);
        laserAttackCD = laserAttackRate;
        rangeAttackCD = rangeAttackRate;
        Invoke("ExitLaserAttack", 1);
    }

    private void ExitLaserAttack()
    {     
        SwitchState(State.ReCharge);
    }

    #endregion

    #region Range Attack

    private void EnterRangeAttackState()
    {
        anim.SetTrigger("isRangeAttack");

    }


    private void BeginRangeAttack()
    {
        Instantiate(armPrefab, rangeAttackPoint.position, Quaternion.identity);
        rangeAttackCD = rangeAttackRate;
        laserAttackCD = laserAttackRate;
        Invoke("ExitRangeAttackState",1.25f);
    }
    private void ExitRangeAttackState()
    {
        SwitchState(State.Idle);
    }

    #endregion


    #region Recharge

    private void EnterReChargeState()
    {
        anim.SetTrigger("isImmune");
        rechargeCD = rechargeRate;
        Invoke("DisableAnimator", 1.95f);
    }

    private void DisableAnimator()
    {
        anim.enabled = false;
    }

    private void UpdateReChargeState()
    {
        if(rechargeCD > 0)
        {
            rechargeCD -= Time.deltaTime;
        }
        else
        {
            ExitReChargeState();
        }
    }

    private void ExitReChargeState()
    {
        anim.enabled = true;
        SwitchState(State.Idle);
    }

    #endregion


    #region Teleport

    private void EnterTeleportState()
    {
        anim.SetTrigger("isTeleport");
    }

    private void UpdateTeleportState()
    {
        if (!isEnraged)
        {
            if (currentPos == 0)
            {
                currentPos = 1;
            }
            else
            {
                currentPos = 0;
            }
            transform.position = positions[currentPos].position;
            SwitchState(State.Idle);
        }
        else
        {
            transform.position = enragePosition.position;
        }
    }
    private void ExitTeleportState()
    {
        if(transform.position.x > FindObjectOfType<PlayerController>().transform.position.x && currentDirection == 1)
        {
            currentDirection = -1;
            transform.Rotate(0, 180, 0f, 0);
        }
        else if(transform.position.x < FindObjectOfType<PlayerController>().transform.position.x && currentDirection == -1)
        {
            currentDirection = 1;
            transform.Rotate(0, 180, 0f, 0);
        }
    }

    #endregion

    #region Enrage State

    private void EnterEnrageState()
    {
        isEnraged = true;
        anim.SetTrigger("isTeleport");
        enrageTime = enrageDuration;
        StartCoroutine(spawnSnowBall());
    }


    private void UpdateEnrageState()
    {
        //Spawn Snow balls for time
       

    }
    private void ExitEnrageState()
    {
        isEnraged = false;
        anim.SetTrigger("isTeleport");
    }
    IEnumerator spawnSnowBall()
    {
        while (enrageTime > 0)
        {
            yield return new WaitForSeconds(snowBallInterval);
            float r = Random.Range(snowBallSpawnPoints[0].position.x, snowBallSpawnPoints[1].position.x);
            Vector2 spawnPoint = new Vector2(r, snowBallSpawnPoints[0].position.y);

            Instantiate(snowBallPrefab, spawnPoint, Quaternion.identity);
            enrageTime -= snowBallInterval;
        }

        SwitchState(State.Idle);
    }

    #endregion



    #region ExtraMethods
    private void SwitchState(State state)
    {
        switch (currentState)
        {
            case State.Idle:                
                break;

            case State.Knockback:
                break;

            case State.MeleeAttack:
                break;
            case State.RangeAttack:
                break;

            case State.LaserAttack:                
                break;

            case State.ReCharge:

                break;
            case State.Teleport:
                ExitTeleportState();
                break;

            case State.Enrage:
                ExitEnrageState();
                break;

            case State.Dead:
                break;

        }

        switch (state)
        {
            case State.Idle:
                EnterIdleState();
                break;

            case State.Knockback:
                break;

            case State.MeleeAttack:
                EnterMeleeAttackState();
                break;

            case State.RangeAttack:
                EnterRangeAttackState();
                break;

            case State.LaserAttack:
                EnterLaserAttackState();
                break;
            case State.ReCharge:
                EnterReChargeState();
                break;

            case State.Teleport:
                EnterTeleportState();
                break;

            case State.Enrage:
                EnterEnrageState();
                break;

            case State.Dead:
                break;

        }
        currentState = state;

    }


    public void Flip()
    {
        facingDirection *= -1;
        transform.Rotate(new Vector3(0, 180.0f, 0));
    }


    #endregion

    private void OnDrawGizmos()
    {

        Gizmos.DrawWireSphere(attackPoint.position, meleeAttackRadius);
        Gizmos.DrawWireSphere(transform.position, rangeAttackDistance);

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
