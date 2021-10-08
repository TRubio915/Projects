using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Tilemaps;

public class MeleeAttackController : MonoBehaviour
{

    private PlayerController controller;
    private RangedAttackWithTeleportation r_AttackController;
    [SerializeField] private GameObject pickUpCollider;
    private Rigidbody2D rb;
    private Animator anim;
    private GameManager manager;

    private float[] attackDetails = new float[2];

    [SerializeField]
    private GameObject bloodParticle;

    [Header("Melee Attack vars")]
    [SerializeField]
    private Transform attackPoint ;
    [SerializeField]
    private Transform upattackPoint;

    [SerializeField]
    private LayerMask whatIsDamageable, whatIsEnvDmgble, whatIsLever;

    public bool
        isUpAttackReady,
        isUpAttacking,
        isAttacking,
        canAttack,
        isAttackingInAir,
        comboAttackBegin;

    private float attackingCoolDown;

    [SerializeField]
    public float
        attackRate,
        attackDamage,
        baseDamage,
        attackRadius,
        attackComboTimer,
        attackComboTimerSet;

    private int
        attackType = -1,
        attackComboCount = 0;

    public float staminaGainOnEnemyHit;

    [SerializeField]
    private GameObject hitParticle;

    [Space(10)]
    [Header("TouchDamage vars")]

    [SerializeField]
    private float touchDamageCD;

    private float
        lastTouchDamageTime;


    [Space(10)]
    [Header("Stop/Slow Time on Attack")]
    [SerializeField]
    private float stopTime;
    [SerializeField]
    private float slowTime;

    public CameraShaker cameraShaker;


    private float _endRange = 0;
    private float _startRange = 255f;
    private float _oscillateRange;
    private float _oscillateOffset;
    private bool bossIsReady = false;

    private int[] values = {0,0,0,0,0,0,0,0,1,1};

    void Start()
    {
        controller = GetComponent<PlayerController>();
        rb = GetComponent<Rigidbody2D>();
        r_AttackController = GetComponent<RangedAttackWithTeleportation>();
        anim = GetComponent<Animator>();
        if (GameObject.FindGameObjectWithTag("Game Manager"))
            manager = GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>();

        isAttacking = false;
        canAttack = true;
        isAttackingInAir = false;
        _oscillateRange = (_endRange - _startRange) / 2;
        _oscillateOffset = _oscillateRange + _startRange;
    }

    // Update is called once per frame
    void Update()
    {

        if (!GameSettings.isGamePaused && !controller.isDead)
        {

            if(Input.GetAxisRaw("Vertical") > 0)
            {
                isUpAttackReady = true;
            }
            else
            {
                isUpAttackReady = false;
            }
            checkIfCanAttack();
            CheckCombatInput();
            UpdateAnimations();
            CheckTouchDamageCoolDown();
        }
    }

    private void CheckCombatInput()
    {
        if (controller.IsPlayerGrounded())
        {
            if (Input.GetButtonDown("Fire1"))
            {
                MeleeAttack();
            }            

        }


    }

    private void checkIfCanAttack()
    {
       
        if (attackingCoolDown <= 0 && r_AttackController.canThrow)
        {
            canAttack = true;
            isAttacking = false;
        }
        else
        {
            attackingCoolDown -= Time.deltaTime;
        }

        if (!isUpAttacking)
        {
            if (comboAttackBegin && attackComboTimer <= 0)
            {
                comboAttackBegin = false;
                attackComboCount = 0;
            }
            else if (attackComboTimer > 0)
            {
                attackComboTimer -= Time.deltaTime;
            }

            attackType = attackComboCount % 3;
        }
        else
        {
            attackComboCount = 0;
            attackType = 3;
        }
    }
    private void MeleeAttack()
    {
        if (canAttack && !r_AttackController.isBeganThrowing)
        {
            if (bossIsReady && Physics2D.OverlapCircle(attackPoint.position, attackRadius, whatIsDamageable))
            {
                int selected = values[Random.Range(0, values.Length)];
                if(selected == 1)
                    FindObjectOfType<WizardBoss>().AttemptToDash();
            }
            r_AttackController.resetThrowingVars();

            if (controller.IsPlayerGrounded())
            {
                rb.velocity = Vector2.zero;
                controller.movementSpeed = 0;
                isAttacking = true;
                if (isUpAttackReady)
                    isUpAttacking = true;
            }
            else
            {
                isAttackingInAir = true;
                //controller.usedUpJump();

            }

            r_AttackController.canThrow = false;
            attackingCoolDown = attackRate;
            attackComboTimer = attackComboTimerSet;
            canAttack = false;
            comboAttackBegin = true;
            attackDamage = baseDamage + 10 * (attackComboCount % 3);
            attackComboCount++;

            FindObjectOfType<AudioManager>().Play("Attack");

            
        }
    }

    private void CheckAttackHitBox()
    {
        

        // Event Called in the animation window for all three attacks
        Collider2D[] detectedObjects;
        if (!isUpAttacking)
            detectedObjects = Physics2D.OverlapCircleAll(attackPoint.position, attackRadius, whatIsDamageable);
        else
            detectedObjects = Physics2D.OverlapCircleAll(upattackPoint.position, attackRadius, whatIsDamageable);

        List<GameObject> objects = new List<GameObject>();
        foreach (Collider2D collider in detectedObjects)
        {
            if (objects.Count > 0)
            {
                if (objects.Contains(collider.gameObject))
                {

                }
                else
                {
                    objects.Add(collider.gameObject);
                }
            }
            else
            {
                objects.Add(collider.gameObject);
            }
        }

        foreach (GameObject obj in objects)
        {
            attackDetails[0] = attackDamage;
            attackDetails[1] = transform.position.x;
            Instantiate(hitParticle, attackPoint.transform.position, Quaternion.Euler(0, 0, Random.Range(0, 360)));
            //collider.transform.SendMessage("SetKnockBackValues", (attackComboCount - 1) % 3);
            obj.transform.SendMessage("Damage", attackDetails);
            manager.IncreasePlayerStamina(staminaGainOnEnemyHit);
            isAttacking = false;
            controller.canMove = true;
            controller.canFlip = true;
            StartCoroutine(cameraShaker.Shake(0.1f, 0.2f));
            FindObjectOfType<AudioManager>().Play("AttackHit");
            //StartCoroutine(stopTimeOnAttack());            
        }
        Collider2D[] detectedLevers = Physics2D.OverlapCircleAll(attackPoint.position, attackRadius, whatIsLever);

        foreach (Collider2D obj in detectedLevers)
        {
            obj.SendMessage("TurnLever");
        }

    }

    IEnumerator stopTimeOnAttack()
    {
        anim.updateMode = AnimatorUpdateMode.Normal;
        GameSettings.isGamePaused = true;
        Time.timeScale = 0;
        yield return new WaitForSecondsRealtime(stopTime);
        GameSettings.isGamePaused = false;
        Time.timeScale = 0.1f;
        yield return new WaitForSecondsRealtime(slowTime);
        Time.timeScale = 1;
        anim.updateMode = AnimatorUpdateMode.UnscaledTime;
    }

    public void Damage(float[] attackDetails)
    {
        if (!controller.isDead && !controller.getUpgradedDashStatus())
        {
            
            int direction;
            int visualDamage = (int)attackDetails[0];

            if (attackDetails[1] < transform.position.x)
            {
                direction = 1;
            }
            else
            {
                direction = -1;
            }

            if (!r_AttackController.isTeleporting)
            {
                GetComponent<RangedAttackWithTeleportation>().resetThrowingVars();
                lastTouchDamageTime = touchDamageCD;
                controller.StopAllBasicMovement();
                controller.KnockBack(direction);
                StartCoroutine(cameraShaker.Shake(0.15f, 0.5f));
                Instantiate(bloodParticle, transform.position, Quaternion.identity);
                FindObjectOfType<AudioManager>().Play("Hurt");
                anim.SetTrigger("isHurt");
                isAttacking = false;
                manager.ReducePlayerHealth(visualDamage);                
                StartCoroutine(stopTimeOnAttack());
                controller.transform.gameObject.layer = 16;
                pickUpCollider.layer = 16;

            }
        }

    }

    private void AttackBegan()
    {
        controller.canMove = false;
        canAttack = false;

    }
    private void AttackFinished()
    {

        controller.canMove = true;
        r_AttackController.canThrow = true;
        canAttack = true;
        isAttacking = false;
        isUpAttacking = false;
        anim.SetBool("isM_attacking", false);
    }

    private void CheckTouchDamageCoolDown()
    {
        if (lastTouchDamageTime > 0)
        {
            lastTouchDamageTime -= Time.deltaTime;
            Color tmp = this.gameObject.GetComponent<SpriteRenderer>().color;

            tmp.a = (_oscillateOffset + Mathf.Sin(Time.deltaTime * 0.35f) * _oscillateRange) / 255;
            this.gameObject.GetComponent<SpriteRenderer>().color = tmp;
        }

        if (lastTouchDamageTime <= 0 && !controller.GetDashStatus())
        {
            this.gameObject.GetComponent<SpriteRenderer>().color = new Color(1, 1, 1, 1);
            controller.transform.gameObject.layer = 8;
            pickUpCollider.layer = 18;
        }
    }

    private void UpdateAnimations()
    {
        anim.SetBool("isM_attacking", isAttacking);
        anim.SetFloat("attackType", attackType);
    }

    private void OnDrawGizmos()
    {
        Gizmos.DrawWireSphere(attackPoint.position, attackRadius);
    }

    private void OnCollisionEnter2D(Collision2D other)
    {
        if (other.collider.CompareTag("Traps"))
        {

            float[] attackDetails = new float[2];
            attackDetails[0] = PlayerStats.PlayerHealth;
            attackDetails[1] = other.collider.transform.position.x;
            if (lastTouchDamageTime <= 0 && !controller.GetDashStatus())
            {
                Damage(attackDetails);
            }
        }
    }

    private void OnCollisionStay2D(Collision2D other)
    {
        if (other.collider.CompareTag("Traps"))
        {

            float[] attackDetails = new float[2];
            attackDetails[0] = 10;
            attackDetails[1] = other.collider.transform.position.x;
            if (lastTouchDamageTime <= 0 && !controller.GetDashStatus())
            {
                Damage(attackDetails);
            }
        }
    }

    public void dmgUpgrade()
    {
        baseDamage += 10;
    }


    public void BossReady()
    {
        bossIsReady = true;
    }
}
